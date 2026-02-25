# FLOOXS Bugs

Six bugs found during FLOOXS evaluation: three in C++ source (`src/BasePDE/`), three in Tcl model procs (`TclLib/Models/Dopant.tcl`). All six are patched in our Docker image.

The C++ bugs were found while investigating the built-in Tcl model procs (`DefectBulk`, `DopantBulk`, `DefectBound`). The Tcl bugs affect the built-in procs when used with the FLOOXS 2026 C++ PDB system.

## C++ bugs (src/BasePDE/)

### 1. Use-after-free on unary + (Reduce.cc)

The expression parser crashes when a PDE equation string begins with `+` (e.g., `" + ddt(X) - ..."`). The crash manifests as heap corruption (`free(): corrupted unsorted chunks`).

`DefectBulk` in `TclLib/Models/Defect.tcl` builds equations by concatenation:

```tcl
proc DefectInit { Mat Sol } {
    pdbSetString $pdbMat $Sol Equation ""
}

proc DefectBulk { Mat Sol } {
    ...
    set eqn "ddt($Sol) - $D0 * grad(Scale$Sol) + IVRecomb"
    set de [pdbGetString $pdbMat $Sol Equation]
    pdbSetString $pdbMat $Sol Equation "$de + $eqn"
}
```

Since `DefectInit` sets the equation to `""`, the concatenation produces `" + ddt(Int) - ..."` with a leading `+`.

**Root cause**: `ExprStore::Unary()` at `Reduce.cc:165`:

```cpp
Expr *ExprStore::Unary( char sign, Expr &ev ) {
    if ( sign == ADD ) return &ev;    // BUG: missing Duplicate()
    // ... other cases all go through MakeMatch (which adds a carrying reference)
}
```

FLOOXS uses reference counting (`Duplicate()` increments, `Deref()` decrements or deletes). Every function returning `Expr*` must return with one "carrying" reference. The ADD branch returns `&ev` without adding one. The caller (`Reduce.cc:747-757`) does `Pop(2)` which derefs the expression (potentially freeing it), then `Push(*e)` which reads through the now-dangling pointer, then `e->Deref()` which writes to freed memory.

**Patch** (`src/BasePDE/Reduce.cc`, line 165):

```cpp
// Before:
if ( sign == ADD ) return &ev;

// After:
if ( sign == ADD ) return ev.Duplicate();
```

**Tcl workaround** (if C++ patch is not applied): check for empty string before concatenation:

```tcl
set de [pdbGetString $pdbMat $Sol Equation]
if {$de eq ""} {
    pdbSetString $pdbMat $Sol Equation $eqn
} else {
    pdbSetString $pdbMat $Sol Equation "$de + $eqn"
}
```

### 2. `term list` requires name (Generic.cc)

`term list` fails with `"Must specify a name to add"`. `DefectBulk` calls `term list` (line 25) to check whether `Pressure{Sol}` already exists.

**Root cause**: In `term_tcl()` (`Generic.cc:249-253`), the name-is-NULL check runs before the list check:

```cpp
char *nm = name->StringValue();
if ( nm == NULL ) {
    tcl.AppendResult("Must specify a name to add", NULL);
    return TCL_ERROR;
}
// ... list check comes later, never reached
```

By contrast, `solution_tcl()` (`Generic.cc:91-97`) correctly handles `list` before the name check.

**Patch** (`src/BasePDE/Generic.cc`, in `term_tcl()`): add a `list` check before the name-is-NULL check:

```cpp
// After the MatSpecified block (~line 248), add:
if ( list->BoolValue() ) {
    int err = Tcl_Eval( interp, "solution list" );
    return err;
}
// Then the existing name check continues
```

### 3. Curly brace stripping corrupts expressions (Generic.cc)

The `solution` command strips all `{` and `}` from stored value strings (`Generic.cc:160-167`):

```cpp
// Eliminate special tcl characters that alagator won't like
if ( v->StringValue() != NULL ) {
  for ( size_t i = 0; i < strlen(v->StringValue()); i++ ) {
    v->StringValue()[i] = (v->StringValue()[i]=='"')?' ':v->StringValue()[i];
    v->StringValue()[i] = (v->StringValue()[i]=='{')?' ':v->StringValue()[i];
    v->StringValue()[i] = (v->StringValue()[i]=='}')?' ':v->StringValue()[i];
  }
}
```

This corrupts two things in expressions stored via `term` (which forwards to `solution`):

1. **Tcl list grouping**: `[DiffLimit Silicon {Int Vac} 0.0]` becomes `[DiffLimit Silicon  Int Vac  0.0]` (4 args instead of 3), causing `wrong # args: should be "DiffLimit Mat Species Barrier"`

2. **Arrhenius prefactor grouping**: `DiffLimit` returns `[Arrhenius {(complex_expr)} act]` where the braces group the first argument. Stripping them breaks the Arrhenius call.

These braces are needed for Tcl argument grouping inside `[...]` expressions. The `Expand()` function in `Parse.cc` evaluates `[...]` via `Tcl_Eval` before the expression parser sees the result, so braces never reach the expression parser ("Alagator").

**Patch** (`src/BasePDE/Generic.cc`, lines 160-167 in `solution_tcl()`): remove the `{` and `}` stripping, keep only `"` stripping:

```cpp
if ( v->StringValue() != NULL ) {
  for ( size_t i = 0; i < strlen(v->StringValue()); i++ ) {
    v->StringValue()[i] = (v->StringValue()[i]=='"')?' ':v->StringValue()[i];
  }
}
```

## Tcl bugs (TclLib/Models/Dopant.tcl)

These three bugs affect the built-in model procs when used with the FLOOXS 2026 C++ PDB system. A patched version of `Dopant.tcl` is in [`../templates/Dopant_patched.tcl`](../templates/Dopant_patched.tcl).

### 4. Segregation key naming (line 206-207)

`Segregation` proc writes interface equations using flat keys:

```tcl
pdbSetString $pdbMat $Sol Equation_$s1 "- $eq"
pdbSetString $pdbMat $Sol Equation_$s2 "$eq"
```

This produces keys like `Oxide_Silicon Boron Equation_Oxide`. But the C++ solver (`BasePDE/Genbc.cc` line 104-138) reads hierarchical paths: `pdbgetString(mn, sn, Mesh1Name, "Equation")` = `Oxide_Silicon Boron Oxide Equation`.

**Fix**: changed to hierarchical keys:

```tcl
pdbSetString $pdbMat $Sol $s1 Equation "- $eq"
pdbSetString $pdbMat $Sol $s2 Equation "$eq"
```

### 5. DopantBulk creates spurious solved variable in oxide (line 31-32)

For `ActiveModel=0` (identity, used for oxide), DopantBulk creates:

```tcl
term name = BoronActive add eqn = Boron Oxide
```

The `term` command maps to `solution name=BoronActive add Oxide solve const val=Boron`. This creates a solved variable `BoronActive` in oxide. The solver tries to include it in the matrix but the gas/oxide boundary has no proper setup for it, producing NaN at Row 1 Col 1.

**Fix**: for `ActiveModel=0`, skip term creation and use `Sol` directly:

```tcl
if {$ActModel == 0} {
    set ActName $Sol       ;# Use "Boron" directly, no term created
} elseif {$ActModel == 1} {
    ...                    ;# Solubility-limited term still created for silicon
}
```

### 6. DopantBulk accesses Charge in neutral materials (line 58)

DopantBulk unconditionally runs `set chg [term name=Charge print $Mat]`, even for `Charge=0` (neutral dopant in oxide). The `term ... print` command maps to `solution name=Charge Oxide print`. If Charge doesn't exist in oxide, the `solution` command creates it (C++ auto-creation in `solution_tcl`, line 134-137 of `BasePDE/Generic.cc`) with no equation, potentially causing solver issues.

**Fix**: guard Charge access behind `chgtype != 0`:

```tcl
if {$chgtype != 0} {
    set chg [term name=Charge print $Mat]
    ...
}
```

## Not fixed: PDB parameter loading

The built-in model procs (`DefectBulk`, `DopantBulk`, etc.) call `pdbDelayDouble Silicon Int Cstar` and similar to look up material parameters. These parameters are defined in the `Params/` directory hierarchy using the old PDB format.

FLOOXS has two PDB implementations:

- **Old** (`Params/pdb` + `Params/pdbHelpers`): Tcl array-based with `__pdbReadParam` for lazy-loading from the `Params/` directory. The parameter files are written for this system. However, it is incomplete: `__pdbGet` (required by `pdbGetDouble` et al.) is defined only in `Params/browser2`, a GUI tool, not in the core helpers.
- **New** (`src/param.tcl`): Thin wrappers around C++ `pdbGet`/`pdbSet` commands (`NameValue.cc`). In-memory tree only. No file-loading mechanism.

`FLOOXS.models` uses the new system (`source $FLXSHOME/src/param.tcl`) and has the old system commented out. As a result, the `Params/` files are never loaded, and any `pdbGetDouble` call for material parameters fails with `"Couldn't find parameter Silicon"`.

Re-enabling the old PDB or writing a loader for the new PDB would require significant work for uncertain benefit. All our simulation templates set parameters and equations explicitly instead of using the built-in procs.

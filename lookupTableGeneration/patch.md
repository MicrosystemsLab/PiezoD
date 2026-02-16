# FLOOXS C++ Patches

Three bugs were found in the FLOOXS C++ source while investigating the
built-in Tcl model procs (`DefectBulk`, `DopantBulk`, `DefectBound`). All
three are in `src/BasePDE/`. These patches are applied but the built-in procs
are still unusable due to the PDB parameter loading issue described at the end
of this document.

---

## Patch 1: Use-after-free on unary `+` (Reduce.cc)

### Problem

The expression parser crashes when a PDE equation string begins with `+`
(e.g., ` + ddt(X) - ...`). The crash manifests as heap corruption:

```
free(): corrupted unsorted chunks
```

### How it happens

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

Since `DefectInit` sets the equation to `""`, the concatenation produces
`" + ddt(Int) - ..."` with a leading `+`.

### Root cause

`ExprStore::Unary()` at `Reduce.cc:165`:

```cpp
Expr *ExprStore::Unary( char sign, Expr &ev ) {
    if ( sign == ADD ) return &ev;    // BUG: missing Duplicate()
    // ... other cases all go through MakeMatch (which adds a carrying reference)
}
```

FLOOXS uses reference counting (`Duplicate()` increments, `Deref()` decrements
or deletes). Every function returning `Expr*` must return with one "carrying"
reference. The ADD branch returns `&ev` without adding one.

The caller (`Reduce.cc:747-757`) does `Pop(2)` which derefs the expression
(potentially freeing it), then `Push(*e)` which reads through the now-dangling
pointer, then `e->Deref()` which writes to freed memory.

### Patch

**File**: `src/BasePDE/Reduce.cc`, line 165

```cpp
// Before:
if ( sign == ADD ) return &ev;

// After:
if ( sign == ADD ) return ev.Duplicate();
```

### Tcl workaround (if C++ patch is not applied)

The crash can be avoided by fixing the Tcl procs that generate the leading `+`:

**`TclLib/Models/Defect.tcl` - `DefectBulk`** (and `DopantDefectPair` in
`Dopant.tcl` line 154, same pattern):

```tcl
# Before (generates leading +):
set de [pdbGetString $pdbMat $Sol Equation]
pdbSetString $pdbMat $Sol Equation "$de + $eqn"

# After (checks for empty):
set de [pdbGetString $pdbMat $Sol Equation]
if {$de eq ""} {
    pdbSetString $pdbMat $Sol Equation $eqn
} else {
    pdbSetString $pdbMat $Sol Equation "$de + $eqn"
}
```

---

## Patch 2 and 3: `term` / `solution` command infrastructure (Generic.cc)

Patches 2 and 3 are in the `term` / `solution` command infrastructure. They
only matter when expressions containing `[DiffLimit ...]` or `[Arrhenius
{...} ...]` are stored via the `term` command, which is what the built-in
model procs do.

### Patch 2: `term list` requires name

#### Problem

`term list` fails with `"Must specify a name to add"`. `DefectBulk` calls
`term list` (line 25) to check whether `Pressure{Sol}` already exists.

#### Root cause

In `term_tcl()` (`Generic.cc:249-253`), the name-is-NULL check runs before
the list check:

```cpp
char *nm = name->StringValue();
if ( nm == NULL ) {
    tcl.AppendResult("Must specify a name to add", NULL);
    return TCL_ERROR;
}
// ... list check comes later, never reached
```

By contrast, `solution_tcl()` (`Generic.cc:91-97`) correctly handles `list`
before the name check.

#### Patch

**File**: `src/BasePDE/Generic.cc`, in `term_tcl()` function

Add a `list` check before the name-is-NULL check (after the `MatSpecified`
block around line 248):

```cpp
// After:
//     Region *m = NULL;
//     if ( pp.MatSpecified() ) m = pp.MaterParse( );
// Add:

    if ( list->BoolValue() ) {
	int err = Tcl_Eval( interp, "solution list" );
	return err;
    }

// Then the existing name check continues:
//     char *nm = name->StringValue();
//     ...
```

### Patch 3: Curly brace stripping corrupts expressions

#### Problem

The `solution` command strips all `{` and `}` from stored value strings
(`Generic.cc:160-167`):

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

This corrupts two things in expressions stored via `term` (which forwards to
`solution`):

1. **Tcl list grouping**: `[DiffLimit Silicon {Int Vac} 0.0]` becomes
   `[DiffLimit Silicon  Int Vac  0.0]` (4 args instead of 3), causing
   `wrong # args: should be "DiffLimit Mat Species Barrier"`

2. **Arrhenius prefactor grouping**: `DiffLimit` returns
   `[Arrhenius {(complex_expr)} act]` where the braces group the first
   argument. Stripping them would break the Arrhenius call.

These braces are needed for Tcl argument grouping inside `[...]` expressions.
The `Expand()` function in `Parse.cc` evaluates `[...]` via `Tcl_Eval` before
the expression parser sees the result, so braces never reach the expression
parser ("alagator").

#### Patch

**File**: `src/BasePDE/Generic.cc`, lines 160-167 in `solution_tcl()`

Remove the `{` and `}` stripping, keep only `"` stripping:

```cpp
// Before:
// Eliminate special tcl characters that alagator won't like
if ( v->StringValue() != NULL ) {
  for ( size_t i = 0; i < strlen(v->StringValue()); i++ ) {
    v->StringValue()[i] = (v->StringValue()[i]=='"')?' ':v->StringValue()[i];
    v->StringValue()[i] = (v->StringValue()[i]=='{')?' ':v->StringValue()[i];
    v->StringValue()[i] = (v->StringValue()[i]=='}')?' ':v->StringValue()[i];
  }
}

// After:
// Eliminate special tcl characters that alagator won't like
// Note: Only strip double quotes. Curly braces must be preserved
// because they are needed for Tcl argument grouping inside [...]
// expressions (e.g. [DiffLimit Silicon {Int Vac} 0.0] and
// [Arrhenius {(complex_expr)} act]). The Expand() function in
// Parse.cc evaluates [...] via Tcl before the expression parser
// sees the result, so braces never reach the expression parser.
if ( v->StringValue() != NULL ) {
  for ( size_t i = 0; i < strlen(v->StringValue()); i++ ) {
    v->StringValue()[i] = (v->StringValue()[i]=='"')?' ':v->StringValue()[i];
  }
}
```

---

## Not fixed: PDB parameter loading

The built-in model procs (`DefectBulk`, `DopantBulk`, etc.) call
`pdbDelayDouble Silicon Int Cstar` and similar to look up material parameters.
These parameters are defined in the `Params/` directory hierarchy (e.g.
`Params/Silicon/Interstitial`, `Params/Silicon/Boron/Info`) using the old PDB
format (`array set $Base {key {Type value}}`).

FLOOXS has two PDB implementations:

- **Old** (`Params/pdb` + `Params/pdbHelpers`): Tcl array-based with
  `__pdbReadParam` for lazy-loading from the `Params/` directory. The parameter
  files are written for this system. However, it is incomplete: `__pdbGet`
  (required by `pdbGetDouble` et al.) is defined only in `Params/browser2`,
  a GUI tool, not in the core helpers.
- **New** (`src/param.tcl`): Thin wrappers around C++ `pdbGet`/`pdbSet`
  commands (`NameValue.cc`). In-memory tree only. No file-loading mechanism.

`FLOOXS.models` uses the new system (`source $FLXSHOME/src/param.tcl`) and
has the old system commented out. As a result, the `Params/` files are never
loaded, and any `pdbGetDouble` call for material parameters fails with
`"Couldn't find parameter Silicon"`.

Re-enabling the old PDB or writing a loader for the new PDB would require
significant work for uncertain benefit. Instead, our simulation scripts set
parameters and equations explicitly (see `test_explicit.tcl`).


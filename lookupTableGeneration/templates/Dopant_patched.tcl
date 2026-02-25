proc InitDopant {Mat Sol} {
    set pdbMat [pdbName $Mat]

    set model [pdbGetSwitch $pdbMat $Sol DiffModel]

    #if we have the react model, we need to set up solutions
    if {$model == 3} {		
	#get a list of defects we react with
	set ld [pdbGetString $pdbMat $Sol Defects]
	foreach d $ld {
	    solution add name = ${Sol}${d} solve !damp !negative
	}
	puts "Reaction Model for $Sol"

    #else we need to make sure things are clear!
    } else {
	#get a list of defects we react with
	set ld [pdbGetString $pdbMat $Sol Defects]
	foreach d $ld {
	    solution name = ${Sol}${d} nosolve !damp !negative
	}
    }
}

proc DopantBulk { Mat Sol } {
    set pdbMat [pdbName $Mat]

    #work with the active model - create DopantActive...
    set ActModel [pdbGetSwitch $pdbMat $Sol ActiveModel]
    set ActName ${Sol}Active
    # FIX: For ActiveModel=0 (identity), don't create a term.
    # "term" maps to "solution ... solve const val=..." which creates a
    # spurious solved variable in the matrix, causing NaN in oxide.
    # Use Sol directly instead.
    if {$ActModel == 0} {
	set ActName $Sol
    } elseif {$ActModel == 1} {
	set ss [pdbDelayDouble $pdbMat $Sol Solubility]
	set eqn "($ss) * $Sol / (($ss) + $Sol)"
	term name = $ActName add eqn = ($eqn) $Mat
    } else {				
	#need to add dynamic precipitation here!
	set ss [pdbDelayDouble $pdbMat $Sol Solubility]
	set eqn "$ss * $Sol / ($ss + $Sol)"
	term name = $ActName add eqn = $eqn $Mat
    }

    set model [pdbGetSwitch $pdbMat $Sol DiffModel]
    set ChgName 0
    if {$model == 0} {		
	set ChgName [DopantConstant $Mat $Sol]
    } elseif {$model == 1} {		
	set ChgName [DopantFermi $Mat $Sol]
    } elseif {$model == 2} {		
	set ChgName [DopantPair $Mat $Sol]
    } elseif {$model == 3} {		
	set ChgName [DopantReact $Mat $Sol]
    }

    #set up charged species in potential equation
    set chgtype [pdbGetSwitch $pdbMat $Sol Charge]
    # FIX: Only access Charge term for charged dopants (chgtype != 0).
    # For neutrals (chgtype=0, e.g. oxide), "term name=Charge print $Mat"
    # creates a spurious Charge solution with no equation, causing NaN.
    if {$chgtype != 0} {
	set chg [term name=Charge print $Mat]
	if {[lsearch $chg $ChgName] == -1} {
	    #acceptor
	    if {$chgtype == 1} {
		term name = Charge add eqn = "$chg - $ChgName" $Mat
	    #donor
	    } elseif {$chgtype == 2} {
		term name = Charge add eqn = "$chg + $ChgName" $Mat
	    }
	}
    }
}

proc DopantConstant { Mat Sol } {
    set pdbMat [pdbName $Mat]
    set diff [pdbDelayDouble $pdbMat $Sol Dstar]
    # FIX: Use Sol directly instead of ${Sol}Active.
    # For ActiveModel=0, DopantBulk sets ActName=Sol (no term created).
    # Using Sol here avoids referencing a non-existent BoronActive variable.
    set eqn "ddt($Sol) - $diff * grad( $Sol )"
    pdbSetString $pdbMat $Sol Equation $eqn
    return $Sol
}


proc DopantFermi { Mat Sol } {
    set pdbMat [pdbName $Mat]

    #buld the diffusivity
    set difnam Diff$Sol
    #build the diffusivity term
    set dif [pdbDelayDouble $Mat $Sol D0]
    if {[pdbIsAvailable $Mat $Sol Dn]} {
	set dif "$dif + [pdbDelayDouble $Mat $Sol Dn] * (Noni)"
    }
    if {[pdbIsAvailable $Mat $Sol Dnn]} {
	set dif "$dif + [pdbDelayDouble $Mat $Sol Dnn] * (Noni)^2"
    }
    if {[pdbIsAvailable $Mat $Sol Dp]} {
	set dif "$dif + [pdbDelayDouble $Mat $Sol Dp] * (Poni)"
    }
    if {[pdbIsAvailable $Mat $Sol Dpp]} {
	set dif "$dif + [pdbDelayDouble $Mat $Sol Dpp] * (Poni)^2"
    }
    term name = $difnam add eqn = $dif $Mat

    set ActName ${Sol}Active
    set eqn "ddt($Sol) - $difnam * grad( $ActName )"
    set chgtype [pdbGetSwitch $pdbMat $Sol Charge]
    if {$chgtype == 1} {
	term name = $difnam add eqn = "($dif) / (Poni)" $Mat
	set eqn "ddt($Sol) - $difnam * grad( $ActName * (Poni) )"
    } elseif {$chgtype == 2} {
	term name = $difnam add eqn = "($dif) / (Noni)" $Mat
	set eqn "ddt($Sol) - $difnam * grad( $ActName * (Noni) )"
    }

    pdbSetString $pdbMat $Sol Equation $eqn
    return $ActName
}


proc DopantDefectPair { Mat Sol Def } {
    puts "DopantDefectPair $Mat $Sol $Def"

    #buld the diffusivity
    set difnam Diff${Sol}${Def}

    #build the diffusivity term
    set dif [pdbDelayDouble $Mat $Sol $Def D0]
    if {[pdbIsAvailable $Mat $Sol $Def Dn]} {
	set dif "$dif + [pdbDelayDouble $Mat $Sol $Def Dn] * Noni"
    }
    if {[pdbIsAvailable $Mat $Sol Dnn]} {
	set dif "$dif + [pdbDelayDouble $Mat $Sol $Def Dnn] * Noni^2"
    }
    if {[pdbIsAvailable $Mat $Sol Dp]} {
	set dif "$dif + [pdbDelayDouble $Mat $Sol $Def Dp] * Poni"
    }
    if {[pdbIsAvailable $Mat $Sol Dpp]} {
	set dif "$dif + [pdbDelayDouble $Mat $Sol $Def Dpp] * Poni^2"
    }
    set SubName ${Sol}Sub
    set chgtype [pdbGetSwitch $Mat $Sol Charge]
    set chg 1.0
    if {$chgtype == 1} {
	set chg Poni
    } elseif {$chgtype == 2} {
	set chg Noni
    }
    puts $dif
    term name = $difnam add eqn = "( $dif ) / $chg" $Mat

    set eqn "$difnam * grad( $SubName * Scale${Def} * $chg )"
    term name = Flux${Sol}${Def} add eqn = $eqn $Mat
    puts "$Sol $eqn"

    set de [pdbGetString $Mat $Sol Equation]
    pdbSetString $Mat $Sol Equation "$de - Flux${Sol}${Def}"
    set de [pdbGetString $Mat $Def Equation]
    set bind [pdbDelayDouble $Mat $Sol $Def Binding]
    pdbSetString $Mat $Def Equation "$de + ddt($bind * $SubName * $Def) - Flux${Sol}${Def}"
}

proc DopantPair { Mat Sol } {
    set pdbMat [pdbName $Mat]

    #get all of the defects we are working with
    set ld [pdbGetString $pdbMat $Sol Defects]

    #build an expression for the substitutional dopant
    set ActName ${Sol}Active
    set den 1.0
    foreach d $ld {
	set den "$den + $d * [pdbDelayDouble $pdbMat $Sol $d Binding]"
    }
    term name = ${Sol}Sub add eqn = "$ActName / ( $den )" $Mat

    #create the basic equation and then add
    pdbSetString $pdbMat $Sol Equation "ddt( $Sol )"

    #for each dopant defect pair, build the flux
    foreach d $ld {
	DopantDefectPair $pdbMat $Sol $d
    }
    puts [pdbGetString $pdbMat $Sol Equation]

    return ${Sol}Sub
}


proc Segregation { Mat Sol } {
    set pdbMat [pdbName $Mat]

    #get the names of the sides
    set s1 [FirstMat $pdbMat]
    set s2 [SecondMat $pdbMat]

    set ss1 [pdbIsAvailable $s1 $Sol DiffModel]
    set ss2 [pdbIsAvailable $s2 $Sol DiffModel]

    if { $ss1 && $ss2 } {
	puts "$pdbMat $Sol"
	set seg [pdbDelayDouble $pdbMat $Sol Segregation]
	set trn [pdbDelayDouble $pdbMat $Sol Transfer]
	puts "$seg $trn"

	set sm1 ${Sol}($s1)
	set sm2 ${Sol}($s2)
	set eq "$trn * ($sm1 - $sm2 / $seg)"
	pdbSetString $pdbMat $Sol $s1 Equation "- $eq"
	pdbSetString $pdbMat $Sol $s2 Equation "$eq"
    }
}



proc DopantDefectReact { Mat Sol Def } {
    puts "DopantDefectReact $Mat $Sol $Def"
    set S ${Sol}${Def}

    #assume the dopant-defect diffusivity is constant
    set B [pdbDelayDouble $Mat $Sol $Def Binding]
    set Cs [pdbDelayDouble $Mat $Def Cstar]
    set D0 [pdbDelayDouble $Mat $Sol $Def D0]
    set dax $D0
    puts $dax

    #do the same charge dependence as for the equilibrium diff
    if {[pdbIsAvailable $Mat $Sol $Def Dn]} {
	set dax "$dax + [pdbDelayDouble $Mat $Sol $Def Dn] * Noni"
    }
    if {[pdbIsAvailable $Mat $Sol Dnn]} {
	set dax "$dax + [pdbDelayDouble $Mat $Sol $Def Dnn] * Noni^2"
    }
    if {[pdbIsAvailable $Mat $Sol Dp]} {
	set dax "$dax + [pdbDelayDouble $Mat $Sol $Def Dp] * Poni"
    }
    if {[pdbIsAvailable $Mat $Sol Dpp]} {
	set dax "$dax + [pdbDelayDouble $Mat $Sol $Def Dpp] * Poni^2"
    }
    set dax "($dax) / ($B * Eq$Def)"
    set SubName ${Sol}Sub
    set chgtype [pdbGetSwitch $Mat $Sol Charge]
    set chg 1.0
    if {$chgtype == 1} {
	set chg Poni
    } elseif {$chgtype == 2} {
	set chg Noni
    }
    puts $dax
    set flux "$dax * grad( $S * $chg ) / $chg"
    puts $flux

    #build the reaction
    set K [pdbDelayDouble $Mat $Sol $Def Krate]
    term name = React$S add eqn = "$K * (${Sol}Sub * $Def - ${Sol}${Def} / $B)" $Mat

    set de [pdbGetString $Mat $Sol Equation]
    pdbSetString $Mat $Sol Equation "$de + React$S"
    set de [pdbGetString $Mat $Def Equation]
    pdbSetString $Mat $Def Equation "$de + React$S"

    pdbSetString $Mat ${Sol}${Def} Equation "ddt($S) - $flux - React$S"
}


proc DopantReact { Mat Sol } {
    set pdbMat [pdbName $Mat]

    #get all of the defects we are working with
    set ld [pdbGetString $pdbMat $Sol Defects]

    #see if we have created a substitutional solution variable

    set ActName ${Sol}Active
    set den $ActName
    foreach d $ld {
	set den "$den - ${Sol}${d}"
    }
    term name = ${Sol}Sub add eqn = "$den" $Mat
    puts "Substitutional is [term name = ${Sol}Sub $Mat print]"

    #create the basic equation and then add to it
    pdbSetString $pdbMat $Sol Equation "ddt( $Sol )"
    puts [pdbGetString $pdbMat $Sol Equation]

    #for each dopant defect pair, build the flux
    foreach d $ld {
	DopantDefectReact $pdbMat $Sol $d
    }
    puts [pdbGetString $pdbMat $Sol Equation]

    return ${Sol}Sub
}

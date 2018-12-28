class binana:

    def _init_
        get receptor filname

        check if ligand was passed as list 
        else filename was passed as ligand 

        assign receptor secondary structure #this is a seperate method 

        VARIABLES TO INIIALIZE by Data Structure: 
            DICTIONARIES: 
                ligand receptor atoms that are less than two and a half angstroms
                ligand receptor atoms that are less than four angstroms 
                electrostatic interactions between receptor and ligand atoms 
                active site flexibility 
                hydrogen bonds
                hydrophobic interactions
                pi electron interactions
                salt bridges 

            PDB: 
                close contacts # if two atoms are < 2.5 A
                contacts # if two atoms are < 4.0 A 
                alpha helix contact # if atoms in contact are < 4.0 A and receptor atom is part of alpha sheet 
                beta sheet contact # if atoms in contact are < 4.0 A and receptor atom is part of beta sheet 
                other 2nd structure contact # if atoms in contact are < 4.0 A and receptor atom is part of neither a beta sheet or alpha helix 
                side chain # if atoms in contact are < 4.0 A and receptor atom is part of a sidechain
                back bone # if atoms in contact are < 4.0 A and receptor atom is part of back bone 
                hydrophobic contacts # if atoms in contact are < 4.0 A and receptor and ligand atoms are both carbons 
                hydrogen bonds # if atoms in contact are < 4.0 A and a hydrogen is at optimal angle to ligand or receptor carbon, it could be h-donor
                pi cation interactions # if atoms in contact are < 6.0 A and are either part of a pi system or are charged then we add to this
                salt bridges # 

        for ligand atom in ligand.AllAtoms dictionary
            for receptor atom in receptor.AllAtoms dictionary
                set ligand atom to a specific atom in all atoms dict # WE NEED TO GO ITERATIVELY THROUGH LIGANDS.ALLATOMS  
                set receptor atom to specific atom in all atom dict
                
                find distance betwen ligand atom and receptor atom 

                if distance is less than 2.5, 
                    then add combinatino of ligand and receptor to "ligand receptor atoms that are less than two and a half angstroms"
                    add receptor and ligand to close contacts
                if distance is less than 4.0, 
                    then add combination of ligand and receptor to "ligand receptor atoms that are less than four angstroms"

                if distance less than 4 
                    calculate electrostatic energies using coulomb energy 
                    add the electrostatic energies to special electrostatic dictionary 

                if distance less than 4 
                    get statistics for active-site flexibility 
                    check if sidechain/backbone 
                    check if atoms are part of beta sheet or alpha helix # BASED ON THESE CRITERIA WE CAN ADD THEM TO DIFFERENT VARIABLES ABOVE

                if distance less than 4
                    look for hydrophobic contacts # LIGAND CARBON AND RECEPTOR CARBON ARE NEXT TO EACH OTHER 
                    add hydrophobic_key to hydrophobics list 
                
                if distance less than 4
                    if ligand or receptor has oxygen and nitrogen 
                        check for hydrogens 
                        check angles of the hydrogens # BASED ON THESE CRITERIA WE ADD TO HYDROGEN BONDS PDB  
        
        get total number of each atom type in ligand 

        for aromatic ring in ligand's aromatic rings 
            for aromatic ring in receptor aromatic rings

                find distance between centers of both aromatic rings

                    if distance less than 7.5 
                        calculate if pi bonds are parallel since this is chemically significant
                        
                        if after math calculatoins, they are roughly parallel
                            re-check probability of pi stacking 
                            
                            for each ligand atom in ligand aromatic ring  
                                project ligand atom onto plane of receptor ring 
                            
                            if pi pi stacking is false 
                                project receptor atom onto plane of ligand ring 
                            
                            if pi pi stacking true 
                                add these structures and atoms to dictionary of pi interactions
                        
                        else if they are roughly perpendicular 
                            iterate through every combo of ligand and receptor atoms to get lowest value found of distance

                            if the lowest distance is <=5 
                                check if ligand pi e's are pointing at receptor pi e's or other way around 
                            
                                if projected point is inside ring of plane being projected on 
                                    update PI (T) interactions PDB
                                
        time to identify pi cation interactions 

        for aromatic ring in receptor 
            for charged atom in ligand charged atoms 
                if positive charge 
                    if charged atom to aromatic center distance less than 6 
                        project charged atom onto aromatic ring plane and save point 
                        if projected charged atom is less than the distance of the atomic radius of the aromatic ring 
                            add this to pi cation interactions PDB 

        for aromatic ring in ligand 
            for charged atom in receptor charged atoms 
                if positive charge 
                    if distance less than 6 A
                        project charged atom onto aromatic ring plane and save point  
                        if projected charged atom is less than the distance of the atomic radius of the aromatic ring 
                            add this to pi cation interactions PDB 

        time to count number of salt bridges

        for receptor charge in receptor charges list 
            for ligand charge in ligand charges list 
                if ligand and receptor have opposite charges 
                    if distance is less than 5.5 
                        add this to salt bridges PDB 

        open vina 

        run vina with given files 

        Print vina output   

        for each found pi interaction event in the pi interactions PDB 
            put stacking, shaped and pi cation interactions into seperate lists 

        create single dictionary and add all the above dictionaries and pdb's to it 

        iterate through every variable of every data structure and check for possible error messages to let the user know of 

        make input vector list and add all the dictionaries that have been created throughout this method into it. # DOES NOT ADD ANY PDB 

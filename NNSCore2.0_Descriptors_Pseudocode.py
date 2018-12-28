{ LIGAND-RECEPTOR ATOMS THAT ARE LESS THAN TWO AND A HALF ANGSTROM'S APART }

Output: {}, maintaining counter of the number of times any combination of atom types are found that are less than 2.5 A away from each other 

for indexable_ligand_atom in ligand.AllAtoms dictionary 
    for indexable_receptor_atom in receptor.AllAtoms dictionary

        find distance between the ligand atom and the receptor atom 
        
        if the distance is less than 2.5 Angstroms
            create [] of ligand atom's atom type and receptor atom's atom type 
            sort this list alphabetically and combine into a string that is the KEY for a hashtable 
            maintain a hashtable {} that we use the previously mentioned KEY to access it's VALUE 
            VALUE = counter that tells the number of times a combination of atom types (KEY) are found that had a distance less than 2.5 A away from each other 

return hashtable {}


{ LIGAND-RECEPTOR ATOMS THAT ARE LESS THAN FOUR ANGSTROM'S APART }

Output: {}, maintaining counter of the number of times any combination of atom types are found that are less than 4 A away from each other 

for indexable_ligand_atom in ligand.AllAtoms dictionary 
    for indexable_receptor_atom in receptor.AllAtoms dictionary

        find distance between the ligand atom and the receptor atom 
        
        if the distance is less than 4 Angstroms
            create [] of ligand atom's atom type and receptor atom's atom type 
            sort this list alphabetically and combine into a string that is the KEY for a hashtable 
            maintain a hashtable {} that we use the previously mentioned KEY to access it's VALUE 
            VALUE = counter that tells the number of times a combination of atom types (KEY) are found that had a distance less than 4 A away from each other 

return hashtable {}


{ ELECTROSTATIC INTERACTIONS BETWEEN RECEPTOR AND LIGAND ATOMS }

Output: {}, maintains total coloumb energy found for any particular ligand/receptor atom type combination (KEY), where the ligand atom and receptor atom are less than 4 A away 

for indexable_ligand_atom in ligand.AllAtoms dictionary 
    for indexable_receptor_atom in receptor.AllAtoms dictionary

        find distance between the ligand atom and the receptor atom 
        
        if the distance is less than 4 Angstroms
            find coloumb energy which is ligand charge * receptor charge / distance * 138.94238460104697e4 # THE LAST FACTOR CONVERTS TO JOULES/MOLE
            create [] of ligand atom type and receptor atom type to be used as KEY (after subsequent steps)
            take previous [], sort alphabetically and join the elements into one string which is the KEY for hashtable
            add this KEY to a hashtable where the VALUE is the total summated coloumb energy of every instance so far of a given KEY 

return hashtable {}


{ ACTIVE SITE FLEXIBILITY }

Output: {}, containing the number of instances of every combination of side-chain/backbone and alpha, beta, and "other" carbons. This is to get statistics on active-site flexibility. 

for indexable_ligand_atom in ligand.AllAtoms dictionary 
    for indexable_receptor_atom in receptor.AllAtoms dictionary

        find distance between the ligand atom and the receptor atom 
        
        if the distance is less than 4 Angstroms
            make a string by finding whether atom is part of "side-chain" or "backbone" add a "_" then find out whether carbon is an "alpha", "beta" or "other" carbon
            take this entire string to use as a KEY to a hashtable 
            add the KEY to the hashtable where the VALUE will be the number of times that this KEY has been encountered 
            PURPOSE of this hashtable is to count number of occurences of any combination of side-chain/backbone and different types of carbons

return hashtable {}


{ HYDROGEN BONDS }

Output: {}, containing number of instances of hydrogen bond donors between ligand and receptor, assorted by KEY. The KEY's are every combination encountered of either "ligand/receptor atom", "Sidechain/Backbone" and "Alpha, Beta and Other" Oxygen and Nitrogen atoms

for indexable_ligand_atom in ligand.AllAtoms dictionary 
    for indexable_receptor_atom in receptor.AllAtoms dictionary

        find distance between the ligand atom and the receptor atom 
        
        if the distance is less than 4 Angstroms

            if ligand atom is Oxygen or Nitrogen AND receptor atom is Oxygen or Nitrogen 

                create [] of hydrogens so we can store all hydrogens encountered that are less than a 1.3 A cutoff from our combination of ligand and receptor oxygen's and nitrogen's

                for atom index in ligand.AllAtoms dictionary
                    if current atom is a hydrogen 
                        if the distance from the current hydrogen and the original ligand atom in the outer for loop < 1.3 A 
                            add a comment to the hydrogen mentioning that it's a ligand hydrogen 
                            append this commented hydrogen to the [] of hydrogens
                            PURPOSE: the purpose of this portion of code is to find any and all potential hydrogen bonds and therefore hydrogen donors

                for atom index in receptor.AllAtoms dictionary
                    if current atom is a hydrogen 
                        if the distance from the current hydrogen and the original receptor atom in the outer for loop < 1.3 A 
                            add a comment to the hydrogen mentioning that it's a receptor hydrogen 
                            append this commented hydrogen to the [] of hydrogens
                            PURPOSE: the purpose of this portion of code is to find any and all potential hydrogen bonds and therefore hydrogen donors

                # we now have a list of possible hydrogen bond donors but we need to refine this list to hydrogens that are more statistically likely to be hydrogen bond donors 

                # the way to refine this list is to check that the angle between the hydrogen, the receptor atom and the ligand atom is less than a certain cutoff

                for a hydrogen in unrefined hydrogens list 
                    if absolute value of (180 - angle between receptor atom, hydrogen and ligand atom * 180.00)/ pi <=40.0 
                        # create KEY to add to hashtable
                        create KEY by creating String that says whether hydrogen is on ligand or receptor, then whether receptor atom is part of Sidechain or Backbone and then receptor atom structure 
                        add this KEY to hydrogen bond donor dictionary, where the dictionary keeps a counter of number of times hydrogen bond donor with specific KEY is encountered

return hashtable {} 


{ HYDROPHOBICS }

Output: {}, containing number of instances receptor atoms of combinations of KEYS (Sidechain/Backbone atoms + different atomic structures) encounter a hydrophobic interaction (where receptor carbon and ligand carbon are less than 4.0 A)  

for indexable_ligand_atom in ligand.AllAtoms dictionary 
    for indexable_receptor_atom in receptor.AllAtoms dictionary

        find distance between the ligand atom and the receptor atom 
        
        if the distance is less than 4 Angstroms
            if ligand atom element is carbon AND receptor atom element is carbon 
                create KEY, which is single string to add to hashtable, where key is whether the receptor atom is part of a Sidechain or Backbone + the receptor atom's structure
                add this key to hydrophobic dictionary which maintains a running counter of number of instances of hydrophobic interactions of each Sidechain/Backbone + receptor atom structures

return hashtable {} 


{ PI ELECTRON INTERACTIONS } 

Output: {}, maintaining counter of pi electron interactions with KEY being every combination of pi-stacking interactions/ pi-T (perpendicular) interactions and various atomic structures encountered 

for ligand_aromatic in ligand.aromatic_rings list 
    for receptor_aromatic in receptor.aromatic_rings list 

        find distance between ligand_aromatic center and receptor_aromatic center 

        if distance is less than 7.5 A # there could be pi pi interactions 

            find angle between the planes # planes being the plane each aromatic ring occupies

            if the planes are close to parallel # it's probably pi pi stacking interaction 

                set boolean flag to say pi pi stacking interaction is false # this is a null hypothesis

                for ligand_ring_index in ligand_aromatic.indices list 

                    project ligand atom onto plane of receptor ring 
                    
                    if the distance from this projection to receptor_aromatic is less than or equal to receptor_aromatic radius + pi padding # pi padding is controversial and used to account for a Van Der Waals radius or for a poor crystal structure 
                        set pi pi flag to true  
                        break 

                if pi pi boolean flag is false # we need to still look at possibilities to conclusively prove the null hypothesis 

                    for receptor_ring_index in receptor_aromatic.indices list 

                        project receptor atom onto plane of ligand ring 

                        if distance of this projection to ligand_aromatic is less than or equal to ligand_aromatic radius + pi padding # pi padding is explained above 
                            set pi pi flag to true 
                            break 
                
                if pi pi flag is true # time to add key and maintain counter in our pi interactions dictionary 

                    find structure of any of the atoms in the receptor_aromatic 
                    if structure is empty string then set structure to be "OTHER" 

                    create KEY to be "STACKING_" + structure

                    add this KEY to pi interaction dictionary/hashtable to maintain counter of every combination encountered of different types of interactions
        
            else if planes are more or less perpendicular # it's probably pi-edge interaction ( pi-T interaction )

                set a min_distance float of 100 to use for comparisons later on 

                    for ligand_index in ligand_aromatic.indices list 
                        set ligand atom to ligand_aromatic.indices[ligand_index]

                        for receptor_index in receptor_aromatic.indices list 
                            set receptor atom to receptor_aromatic.indices[receptor_index] 
                            
                            find distance between ligand and receptor atom 
                            
                            if found distance is less than min_distance 
                                set min_distance to be found distance # we do this to keep finding the lowest possible distance encoutnered over all combinations of ligand and receptor atoms 
                
                if min_distance < = 5 # the rings come within less than 5 A at their closest points 

                    # we need to find out which way the pi electrons of either receptor or ligand are pointing 
                    
                    project ligand_aromatic center onto receptor_aromatic
                    project receptor_aromatic center onto ligand_aromatic

                    # a pi-T interaction would mean that projected point should fall within the ring it's projected onto 

                    if either of the points actually fall within the ring onto which they were originially projected onto 

                        find structure of any of the atoms in the receptor_aromatic 
                        if structure is empty string then set structure to be "OTHER" 

                        create KEY to be "T-SHAPED_" + structure

                        add this KEY to pi interaction dictionary/hashtable to maintain counter of every combination encountered of different types of interactions

return hashtable {}


{ SALT BRIDGES }

Output: {}, maintaining counter of number of times salt bridge is encountered for every atomic structure type 

for receptor_charge in receptor.charges list 
    for ligand_charge in ligand.charges list 
        if ligand and receptor charges have oppposite charges # this is the only way you will encounter a salt-bridge complex possibility  
            if distance between charged ligand and receptor atom is less than 5.5 
                
                get structure of receptor atom that is charged 
                if the structure is empty string then set structure to be "OTHER"

                create KEY to be "SALT-BRIDGE_" + structure 

                add key to salt bridges dictionary/hashtable to maintain counter of instances of salt bridge occurence for each atomic structure type 

return hashtable {}
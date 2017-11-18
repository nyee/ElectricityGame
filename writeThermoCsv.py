from rmgpy import settings
import pandas as pd
import re
import os
from rmgpy.species import Species
from rmgpy.data.rmg import RMGDatabase
from rmgpy.chemkin import saveSpeciesDictionary, loadSpeciesDictionary, getSpeciesIdentifier

def addGAVComparison(validation_test_df, H298s_gav, Cp300s_gav, S298s_gav, comments_gav):

    #compare Cp
    validation_test_df['Cp300_gav(cal/mol/K)'] = pd.Series(Cp300s_gav, index=validation_test_df.index)
    gav_lib_diff = abs(validation_test_df['Cp300_gav(cal/mol/K)']-validation_test_df['Cp300_lib(cal/mol/K)'])
    validation_test_df['Cp300_gav_lib_diff(cal/mol/K)'] = pd.Series(gav_lib_diff, index=validation_test_df.index)

    #repeat for Hf
    validation_test_df['H298_gav(kcal/mol/K)'] = pd.Series(H298s_gav, index=validation_test_df.index)
    gav_lib_diff1 = abs(validation_test_df['H298_gav(kcal/mol/K)']-validation_test_df['H298_lib(kcal/mol/K)'])
    validation_test_df['H298_gav_lib_diff(kcal/mol/K)'] = pd.Series(gav_lib_diff1, index=validation_test_df.index)

    #repeat for Sf
    validation_test_df['S298_gav(cal/mol/K)'] = pd.Series(S298s_gav, index=validation_test_df.index)
    gav_lib_diff2 = abs(validation_test_df['S298_gav(cal/mol/K)']-validation_test_df['S298_lib(cal/mol/K)'])
    validation_test_df['S298_gav_lib_diff(cal/mol/K)'] = pd.Series(gav_lib_diff2, index=validation_test_df.index)

    #save comments
    validation_test_df['GAV comments'] = pd.Series(comments_gav, index = validation_test_df.index )

    return validation_test_df
#######################################################
if __name__ == "__main__":

    #get all library names
    libraries=os.listdir(os.path.join(settings['database.directory'], 'thermo/libraries'))
    libraries = [re.sub('\.py', '', name) for name in libraries if re.search('\.py', name)]
    print libraries
    database = RMGDatabase()

    database.load(settings['database.directory'], thermoLibraries=libraries,\
                 kineticsFamilies='none', kineticsDepositories='none', reactionLibraries = [])

    thermoDatabase = database.thermo

    validation_test_dict = {} # key: spec.label, value: (thermo_gav, thermo_lib)

    #filter out some of the species that come from libraries
    bannedSmilesRegex = '[lNSAe45\+\-]' #filter for N,S, +, -, A for Ar, l for Cl, and e for He/Ne


    spec_labels = []
    spec_dictionary = {}
    H298s_lib = []
    Cp300s_lib = []
    S298s_lib = []
    H298s_gav = []
    Cp300s_gav = []
    S298s_gav = []
    comments_gav = []
    count = 0
    for libraryname, library in thermoDatabase.libraries.iteritems():
        for entryName, entry in library.entries.iteritems():
            if entry.data is None:
                print libraryname, entryName
                continue
            smiles = entry.item.toSMILES()
            if re.search(bannedSmilesRegex, smiles): continue

            #save to spec_dictionary
            species = Species()
            species.molecule.append(entry.item)
            species.generateResonanceIsomers()
            baseName = getSpeciesIdentifier(species)
            newName = baseName
            while newName in spec_dictionary.keys():
                newName = baseName + "{0}".format(count)
                count+=1
            if len(newName) > 15: newName = "Renamed" + str(count)
            species.index = -1
            species.label = newName
            spec_dictionary[newName] = species

            spec_labels.append(newName)
            H298s_lib.append(entry.data.getEnthalpy(298)/4184.0) # unit: kcal/mol-K
            S298s_lib.append(entry.data.getEntropy(298)/4.184) # unit: cal/mol-K
            Cp300s_lib.append(entry.data.getHeatCapacity(300)/4.184) # unit: cal/mol-K


            thermo_gav = thermoDatabase.getThermoDataFromGroups(species)
            H298s_gav.append(thermo_gav.getEnthalpy(298)/4184.0) # unit: kcal/mol-K
            Cp300s_gav.append(thermo_gav.getHeatCapacity(300)/4.184) # unit: cal/mol-K
            S298s_gav.append(thermo_gav.getEntropy(298)/4.184) #  unit: cal/mol-K
            comments_gav.append(thermo_gav.comment)

    new_validation_test_df = pd.DataFrame(index=spec_labels)
    new_validation_test_df['Cp300_lib(cal/mol/K)'] = pd.Series(Cp300s_lib, index=new_validation_test_df.index)
    new_validation_test_df['H298_lib(kcal/mol/K)'] = pd.Series(H298s_lib, index=new_validation_test_df.index)
    new_validation_test_df['S298_lib(cal/mol/K)'] = pd.Series(S298s_lib, index=new_validation_test_df.index)

    new_validation_test_df = addGAVComparison(new_validation_test_df, H298s_gav, Cp300s_gav, S298s_gav, comments_gav)

    with open('/Users/Nate/Dropbox (MIT)/Research/RMG/thermo/comparison_scripts/new.csv', 'w') as fout:
        new_validation_test_df.to_csv(fout)
    saveSpeciesDictionary('/Users/Nate/Dropbox (MIT)/Research/RMG/thermo/comparison_scripts/species_dictionary.txt', spec_dictionary.values())

    test = loadSpeciesDictionary('/Users/Nate/Dropbox (MIT)/Research/RMG/thermo/comparison_scripts/species_dictionary.txt')
    for label in spec_labels:
        assert label in test, label
    ###################################################################################################
    # #now get GAVS for benchmark
    # H298s_gav = []
    # Cp300s_gav = []
    # S298s_gav = []
    # comments_gav = []
    #
    #
    # #reload database
    # database = RMGDatabase()
    #
    # database.load(settings['database.directory'], thermoLibraries=libraries,\
    #              kineticsFamilies='none', kineticsDepositories='none', reactionLibraries = [])
    # thermoDatabase = database.thermo
    #
    # for smiles in spec_labels:
    #     species = Species().fromSMILES(smiles)
    #     thermo_gav = thermoDatabase.getThermoDataFromGroups(species)
    #     H298s_gav.append(thermo_gav.getEnthalpy(298)/4184.0) # unit: kcal/mol-K
    #     Cp300s_gav.append(thermo_gav.getHeatCapacity(298)/4.184) # unit: cal/mol-K
    #     S298s_gav.append(thermo_gav.getEntropy(298)/4.184) #  unit: cal/mol-K
    #     comments_gav.append(thermo_gav.comment)
    #
    # old_validation_test_df = addGAVComparison(old_validation_test_df, H298s_gav, Cp300s_gav, S298s_gav, comments_gav)
    #
    # with open('/Users/Nate/Dropbox (MIT)/Research/RMG/thermo/comparison_scripts/old.csv', 'w') as fout:
    #     old_validation_test_df.to_csv(fout)

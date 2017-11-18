from rmgpy import settings
import pandas as pd
from rmgpy.species import Species
from rmgpy.data.rmg import RMGDatabase
import csv
from rmgpy.chemkin import loadSpeciesDictionary

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

if __name__ == "__main__":
    #read csv file and take all library data
    inputCsv = "/Users/Nate/Dropbox (MIT)/Research/RMG/thermo/comparison_scripts/new.csv"
    dictionaryPath = "/Users/Nate/Dropbox (MIT)/Research/RMG/thermo/comparison_scripts/species_dictionary.txt"
    spec_dictionary = loadSpeciesDictionary(dictionaryPath)

    specColumn = 0
    cpColumn = 1
    hColumn = 2
    sColumn = 3

    spec_labels = []
    H298s_lib = []
    Cp300s_lib = []
    S298s_lib = []
    with open(inputCsv) as csvfile:
        spamreader = csv.reader(csvfile)
        for index, row in enumerate(spamreader):
            if index == 0: continue
            spec_labels.append(row[specColumn])
            Cp300s_lib.append(float(row[cpColumn]))
            H298s_lib.append(float(row[hColumn]))
            S298s_lib.append(float(row[sColumn]))

    old_validation_test_df = pd.DataFrame(index=spec_labels)
    old_validation_test_df['Cp300_lib(cal/mol/K)'] = pd.Series(Cp300s_lib, index=old_validation_test_df.index)
    old_validation_test_df['H298_lib(kcal/mol/K)'] = pd.Series(H298s_lib, index=old_validation_test_df.index)
    old_validation_test_df['S298_lib(cal/mol/K)'] = pd.Series(S298s_lib, index=old_validation_test_df.index)

    #now get GAVS for benchmark
    H298s_gav = []
    Cp300s_gav = []
    S298s_gav = []
    comments_gav = []

    #reload database
    database = RMGDatabase()

    database.load(settings['database.directory'], thermoLibraries=[],\
                 kineticsFamilies='none', kineticsDepositories='none', reactionLibraries = [])
    thermoDatabase = database.thermo

    for label in spec_labels:
        species = spec_dictionary[label]
        species.generateResonanceIsomers()
        thermo_gav = thermoDatabase.getThermoDataFromGroups(species)
        H298s_gav.append(thermo_gav.getEnthalpy(298)/4184.0) # unit: kcal/mol-K
        Cp300s_gav.append(thermo_gav.getHeatCapacity(298)/4.184) # unit: cal/mol-K
        S298s_gav.append(thermo_gav.getEntropy(298)/4.184) #  unit: cal/mol-K
        comments_gav.append(thermo_gav.comment)

    old_validation_test_df = addGAVComparison(old_validation_test_df, H298s_gav, Cp300s_gav, S298s_gav, comments_gav)

    with open('/Users/Nate/Dropbox (MIT)/Research/RMG/thermo/comparison_scripts/old.csv', 'w') as fout:
        old_validation_test_df.to_csv(fout)
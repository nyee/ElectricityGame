import pylab
import numpy
import csv
import sys
from IPython.display import display

def clean_headers(headers):
  cleaned_headers = []
  del headers[0]  # remove time header
  for line in headers:
    tokens = line.replace('>','').replace('<','').split(' ')
    del tokens[0]; del tokens[0]; del tokens[-1]
    cleaned_headers.append(': '.join(tokens))
  return cleaned_headers

top_reactions = []
if __name__ == "__main__":

    filename = "/Users/Nate/Dropbox (MIT)/Research/Peng/cresols/merge_v2/merge_v204_OH_sensitivity.csv"
    outputFile = "/Users/Nate/Dropbox (MIT)/Research/Peng/cresols/merge_v2/sensitivity.png"
    # Choose temperature as first command line argument

    targetTime = 0.036

    # for temp in range(len(temperatures)):
      # Choose run to analyze
    # temperature = temperatures[temp]	# temperature in K

    # Extract headers and data
    f = file(filename)
    reader = csv.reader(f)
    headers = reader.next()
    f.close()
    headers = clean_headers(headers)
    data = numpy.loadtxt(open(filename,"rb"),delimiter=",",skiprows=1)
    time = data[:,0]


    #find time row
    for index, timepoint in enumerate(time):
        if timepoint > targetTime:
            row = index -1
            break



    values = []
    for i in range(len(headers)):
      sensitivity = 0
      for j in range(len(time)-1):
        dt = time[j+1]-time[j]
        sensitivity += data[j,i+1]*dt
      values.append((headers[i],sensitivity))

    # Sort sensitivities in descending order
    values = sorted(values, key=lambda value: abs(value[1]))

    reactions = zip(*values)[0][-10:]
    sensitivities = zip(*values)[1][-10:]



    fig = pylab.figure(figsize=(12,8))
    #subplot_num= '22' + str(temp)
    #pylab.subplot(int(subplot_num))
    position = numpy.arange(len(reactions)) + .5  	# bar centers on y-axis
    pylab.barh(position, sensitivities, align='center')
    pylab.yticks(position, reactions, fontsize=8)
    pylab.gcf().subplots_adjust(left=0.5)
    pylab.xlabel('Sensitivity')
    # pylab.title('DIPK OH Sensitivity at ' + str(temperature) + 'K')
    # figure_name = 'figures/' + str(temperature) + 'K.png'
    pylab.title('DIPK OH Sensitivity')
    figure_name = 'figures/' + 'sensitivity.png'
    fig.show()
    # fig.savefig(figure_name)

    reactions = list(reactions)
    reactions.reverse()
    for rxn in reactions:
      if not rxn in top_reactions:
        top_reactions.append(rxn)

      print 'Most Sensitive Reactions'
      for rxn in top_reactions:
        print rxn


(forked from tpmccauley/dimuon-filter)

This package is a filter that selects events from the MultiJet primary dataset from the CMS open
data release, and computes the razor variables MR and Rsq, used in
supersymmetric particle searches. In particular, the event is selected
if there are at least two "Paticle Flow" jets with pT>30 GeV and |eta|<3. A csv
file and a ROOT file containing the MR, Rsq, HT, MET, number of jets, and number of b-jets for each event is produced.

See http://opendata.cern.ch for more information and for context on the instructions below.

To produce files in the VM open a terminal with the X terminal emulator (an icon bottom-left of the VM screen)
and input the commands as explained below.

* Create a CMSSW environment: 

```
    cmsrel CMSSW_4_2_8
```

* Change to the CMSSW_4_2_8/src/ directory:

```
    cd CMSSW_4_2_8/src/
```
* Initialize the CMSSW environment:

```
    cmsenv
```
* Clone the source code:

```
    git clone https://github.com/jmduarte/SUSYBSMAnalysis-RazorFilter SUSYBSMAnalysis/RazorFilter
````
* Compile the code with the command:

```
    scram b
```
* Go to the source directory:

```
    cd SUSYBSMAnalysis/RazorFilter
```
* Run the example configuration file (see comments in the file on changing parameters):

```
    cmsRun MultiJetRun2010B.py
```
which will produce a csv file: MultiJetRun2010B.csv and a ROOT file:
    MultiJetRun2010B.root. Enjoy!

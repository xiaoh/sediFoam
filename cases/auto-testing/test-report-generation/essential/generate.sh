# set git repository directory
codeDir="/l/sunrui/GIT/solvers/lammpsfoamfine/lammpsfoamFineCode"

# No need to edit below
rm lammps-paper.tex 
touch lammps-paper.tex

reportDir=$PWD
cd $codeDir
codeHash=`git log -1 --format="%h"`
cd $reportDir

cat body/before.tex > lammps-paper.tex
versionString="The code version is: "
versionString+=$codeHash
echo $versionString >> lammps-paper.tex

./singleParticle.py
./multiParticle.py
./pressureDrop.py
./expMueller.py

cat body/after.tex >> lammps-paper.tex

make
version=`uname`
# Use different options according to different versions
if [ $version == "Linux" ]
then
    echo "The version you choose is openmpi version"
    evince lammps-paper.pdf
elif [ $version == "Darwin" ]
then
    echo "The version you choose is mac version"
    open lammps-paper.pdf
else
    echo "Sorry, we haven't got the required version."
    echo "Please contact the developer (sunrui@vt.edu) for help."
fi

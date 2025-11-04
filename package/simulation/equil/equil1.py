import os
import importlib

utils = importlib.import_module('..utils', 'utilities.')


class equil1:
    def __init__(self, vars):
        self.__dict__ = vars

    def run(self):
        print('''\n——Running equil1''')
        if not os.path.exists('equil1'):
            os.mkdir('equil1')
        os.chdir('equil1')
        
        try:
            check = self.check_trj_len.check(self.solvprmtop, 'equil1.dcd', self.eq1len)
        except Exception:
            check = False

        # ADDED: get cell vectors and origin from VMD ##########
        if not os.path.exists('equil1.dcd') or check == False:
            with open("get_celldimension.vmd", 'w') as f:
                f.write(f"""mol delete all;
        mol load parm7 {self.solvprmtop} pdb {self.solvpdb}
        set all [atomselect top all];
        set origin [measure center $all];
        set box [measure minmax $all];
        set min [lindex $box 0];
        set max [lindex $box 1];
        set cell [vecsub $max $min];
        set cellx [lindex $cell 0]
        set celly [lindex $cell 1]
        set cellz [lindex $cell 2]
        set cellBasisVector1 [list $cellx 0.0 0.0]
        set cellBasisVector2 [list 0.0 $celly 0.0]
        set cellBasisVector3 [list 0.0 0.0 $cellz]
        put "center, $origin"
        put "cellBasisVector1, $cellBasisVector1"
        put "cellBasisVector2, $cellBasisVector2"
        put "cellBasisVector3, $cellBasisVector3"
        quit""")

            os.system(f"{self.vmd_path} -dispdev text -e get_celldimension.vmd > celldimension.log")

            with open("celldimension.log",'r') as f:
                lines = f.readlines()
                for line in lines:
                    if line.startswith('center'):
                        origin = line.split(',')[1].rstrip('\n').lstrip(' ')
                for line in lines:
                    if line.startswith('cellBasisVector1'):
                        vector1 = line.split(',')[1].rstrip('\n').lstrip(' ')
                for line in lines:
                    if line.startswith('cellBasisVector2'):
                        vector2 = line.split(',')[1].rstrip('\n').lstrip(' ')
                for line in lines:
                    if line.startswith('cellBasisVector3'):
                        vector3 = line.split(',')[1].rstrip('\n').lstrip(' ')
            
            with open("equil1.nvt", 'w') as f:
                f.write(f"""# configuration file for equil1
amber              on
parmfile           {self.solvprmtop}
coordinates        {self.solvpdb}
exclude            scaled1-4
oneFourScaling     0.833333
scnb               2
switching          on
switchdist         10.5
cutoff             12.0
pairlistdist       14.0
temperature        {self.T_start}
cellBasisVector1   {vector1}
cellBasisVector2   {vector2} 
cellBasisVector3   {vector3}
cellOrigin         {origin}
outputName         equil1
dcdUnitCell        yes
dcdFreq            {self.dcdfreq}
restartfreq        500
outputEnergies     500
outputPressure     500
outputtiming       500
XSTFreq            500
binaryoutput       no
binaryrestart      no
hgroupcutoff       2.8
wrapAll            off
wrapWater          on
langevin           on
langevinTemp       {self.T_start}
langevinDamping    1
langevinHydrogen   no
PME                yes
PMEGridSpacing     1.0
PMETolerance       10e-6
PMEInterpOrder     4
timestep           {self.timestep}
fullelectfrequency 2
nonbondedfreq      1
rigidbonds         all
rigidtolerance     0.00001
rigiditerations    400
stepspercycle      10
splitpatch         hydrogen
margin             2
useflexiblecell    no
useConstantRatio   no
fixedAtoms         on
fixedAtomsForces   on 
fixedAtomsFile     {self.solvrest}
fixedAtomsCol      B
minimize           {self.minsteps}
reinitvels         {self.T_start}
run                {self.eq1len}

            os.system("rm get_celldimension.vmd celldimension.log")
            os.system(f"{self.engine} +p{self.n_procs} +setcpuaffinity +devices {self.device} equil1.conf | tee equil1.log | grep 'ENERGY'")
        ########################################################

        eq1out = {
            'coor': os.path.abspath('output.coor'),
            'vel': os.path.abspath('output.vel'),
            'xsc': os.path.abspath('output.xsc')
            }

        update = {'eq1': eq1out}

        self.output = {}
        self.output |= update

        os.chdir('..')

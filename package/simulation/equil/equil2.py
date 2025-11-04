import os
import MDAnalysis as mda
import importlib



class equil2:
    def __init__(self, vars):
        self.__dict__ = vars
        
    def run(self):
        print('\n——Running equil2')
        if not os.path.exists('equil2'):
            os.mkdir('equil2')
        os.chdir('equil2')
        
        if os.path.exists('equil2.dcd'):
            try:
                check = self.check_trj_len.check(self.solvprmtop, 'equil2.dcd', self.eq2len)
            except Exception:
                check = False

        if not os.path.exists('equil2.dcd') or check == False:
            out = self.output['eq1']

        # ADDED: configuration for NAMD ######
            with open("equil2.npt", 'w') as f:
                f.write(f"""# configuration file for equil2
amber                  on
parmfile               {self.solvprmtop}
coordinates            {out['coor']}
velocities             {out['vel']}
extendedSystem         {out['xsc']}
exclude                scaled1-4
oneFourScaling         0.833333
scnb                   2
switching              on
switchdist             10.5
cutoff                 12.0
pairlistdist           14.0
outputName             equil2
dcdUnitCell            yes
dcdFreq                {self.dcdfreq}
restartfreq            500
outputEnergies         500
outputPressure         500
outputtiming           500
XSTFreq                500
binaryoutput           yes
binaryrestart          yes
hgroupcutoff           2.8
wrapAll                off
wrapWater              on
langevin               on
langevinTemp           {self.T_start}
langevinDamping        1
langevinHydrogen       no
langevinPiston         on
langevinPistonTarget   1.01325
langevinPistonPeriod   200
langevinPistonDecay    100
langevinPistonTemp     {self.T_start}
useflexiblecell        yes
useConstantRatio       yes
ExcludeFromPressureFile {self.solvrest}
ExcludeFromPressureCol B
useGroupPressure       yes
PME                    yes
PMEGridSpacing         1.0
PMETolerance           10e-6
PMEInterpOrder         4
timestep               {self.timestep}
fullelectfrequency     2
nonbondedfreq          1
rigidbonds             all
rigidtolerance         0.00001
rigiditerations        400
stepspercycle          10
splitpatch             hydrogen
margin                 2
fixedAtoms             on
fixedAtomsForces       on 
fixedAtomsFile         {self.solvrest}
fixedAtomsCol          B
run                    {self.eq2len}
""")

            os.system(f"{self.engine} +p{self.n_procs} +setcpuaffinity +devices {self.device} equil2.conf | tee equil2.log | grep 'ENERGY'")
        #################################
        
        if os.path.exists('wrap.dcd'):
            try:
                check = self.check_trj_len.check(self.solvprmtop, 'wrap.dcd', self.eq2len)
            except Exception:
                check = False

        if not os.path.exists('wrap.dcd') or check == False:
            self.wrapping.wrap_equil2(self.solvpdb, 'equil2.dcd', 'wrap.dcd')

        if not os.path.exists('eq2_last.pdb'):
            wrap_u = mda.Universe(self.solvpdb, 'wrap.dcd')
            wrap_u.trajectory[-1]

            with mda.Writer('eq2_last.pdb', wrap_u.atoms.n_atoms) as W:
                W.write(wrap_u.atoms)

        eq2out = {
            'dcd': os.path.abspath('wrap.dcd'),
            'coor': os.path.abspath('output.coor'),
            'vel': os.path.abspath('output.vel'),
            'xsc': os.path.abspath('output.xsc')
            }

        update = {'eq2': eq2out}

        solvpdb = os.path.abspath('eq2_last.pdb')
        self.solvpdb = solvpdb

        self.output|= update

        os.chdir('..')

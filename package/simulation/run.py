import os
import sys
import importlib



class simulation:
    def __init__(self, vars):
        self.__dict__ = vars



    def run(self):
        print(f'\n——Running simulation\n')
        if not os.path.exists('MD'):
            os.mkdir('MD')

        os.chdir('MD')

        self.done_temp = []

        trj_list = []
        simulation_length = 0

        for i,set in enumerate(self.temperature):
            simulation_length += set[1]

            trj = self.run_simulation(i,set)
            wrap_trj = self.wrap_simulation(trj, set)

            if self.dryer == 'yes':
                dry_trj = self.wrapping.dry_trj(wrap_trj, set)
                traj = os.path.abspath(dry_trj)

            elif self.dryer == 'no':
                traj = os.path.abspath(wrap_trj)

            avg = self.score_simulation(trj, i, set)

            trj_list.append(traj)

            self.done_temp.append(set[0])

            if avg > self.score_stop:
                print(f'    {avg} > {self.score_stop}')
                print('    STOP')
                break

                self.t_end = set[0]

            else:
                continue

        self.final_dcd = self.merge_run(trj_list, simulation_length)

        os.chdir('..')

        



    def run_simulation(self, i, set):
        temp, length = set

        if not os.path.exists('output_files'):
            os.mkdir('output_files')

        if not os.path.exists(f'run_{temp}.dcd'):
            print(f'    run_{temp}.dcd doesn\'t exists')

            self.run_temp(i, temp, length)

        else:
            print(f'    run_{temp}.dcd found')

            lgt = int(length) / int(self.dcdfreq)
            check = self.check_trj_len.check(self.solvprmtop, f'run_{temp}.dcd', gt)
            
            if check == False:
                if self.resume == True:
                    print('    Resuming trajectory')
                else:
                    print('    Restarting trajectory')

                self.run_temp(i, temp, length)

        xsc = os.path.abspath(f'output_files/output_{temp}.xsc')
        coor = os.path.abspath(f'output_files/output_{temp}.coor')
        vel = os.path.abspath(f'output_files/output_{temp}.vel')
        
        output = {
                'xsc': xsc,
                'coor': coor,
                'vel': vel
                }

        self.output[i] = output

        return f'run_{temp}.dcd'


    
    def run_temp(self, i, temp, length):
        # #### first step starts from pdb (last frame of equil2), prmtop and xsc, as in sumd 

        if i == 0:
            dict = 'eq2'
        else:
            dict = i - 1

        xsc = self.output[dict]['xsc']
        coor = self.output[dict]['coor']
        vel = self.output[dict]['vel']

        # ADDED: configuration needed to work with NAMD ##
        if i == 0:
            input_files = f'''extendedSystem     {xsc}
temperature        {temp}'''

        else:
            input_files = f'''extendedSystem     {xsc}
binCoordinates     {coor}
binVelocities      {vel}'''

        with open("run.conf", 'w') as f:
                f.write(f"""# configuration for main run
amber              on
parmfile           {self.solvprmtop}
coordinates        {self.solvpdb}
{input_files}
exclude            scaled1-4
oneFourScaling     0.833333
scnb               2
switching          on
switchdist         10.5
cutoff             12.0
pairlistdist       14.0
dcdUnitCell        yes
dcdFreq            {self.dcdfreq}
restartfreq        500
outputEnergies     500
outputPressure     500
outputtiming       500
XSTFreq            500
binaryoutput       yes
binaryrestart      yes
hgroupcutoff       2.8
wrapAll            off
wrapWater          on
langevin           on
langevinTemp       {temp}
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
outputname         output
dcdfile            run_{temp}.dcd
run                {length}
""")

        os.system(f"{self.engine} +p{self.n_procs} +setcpuaffinity +devices {self.device} run.conf | tee run.log | grep 'ENERGY'")
        ##################################################

        os.system(f"cp output.coor output_files/output_{temp}.coor")
        os.system(f"cp output.vel output_files/output_{temp}.vel")
        os.system(f"cp output.xsc output_files/output_{temp}.xsc")
        os.system('rm -r restart*')



    def wrap_simulation(self, trj, set):
        temp, length = set    
        wrap_trj = f"swag_{temp}.dcd"
        
        if not os.path.exists(wrap_trj):
            print(f'\n    {wrap_trj} doesn\'t exists')
            self.wrapping.run(self.solvpdb, trj, wrap_trj)

        else:
            lgt = int(length) / int(self.dcdfreq)
            check = self.check_trj_len.check(self.solvprmtop, wrap_trj, lgt)
            if check == False:
                self.wrapping.run(self.solvprmtop, self.solvpdb, trj, wrap_trj)

        wt = os.path.abspath(wrap_trj)

        return wt



    def score_simulation(self, wrap_trj, i, set):
        temp, length = set

        if not os.path.exists(f'score_{temp}'):
            outscore = self.score(self.complprmtop, wrap_trj, temp)

            with open(f'score_{temp}', 'w') as f:
                for s in outscore:
                    f.write(f'{s}\n')

        else:
            outscore = []

            with open(f'score_{temp}', 'r') as f:
                lines = f.readlines()
                for line in lines:
                    outscore.append(float(line.rstrip('\n')))

        n = int(len(outscore)*self.stop_range/100)
        stoplist = outscore[-n:]
        
        stop_avg = calc_avg(stoplist)

        output = {'results': outscore,
                'avg': stop_avg,
                }

        self.output[i] |= output

        return stop_avg



    def merge_run(self, trj_list, length):
        finaltrj = 'simulation.dcd'

        if self.dryer == 'yes':
            topology = self.complprmtop
        elif self.dryer == 'no':
            topology = self.solvprmtop

        if not os.path.exists(finaltrj):
            check = False

        else:
            lgt = int(length) / int(self.dcdfreq)
            check = self.check_trj_len.check(topology, finaltrj, lgt)
        
        self.wrapping.merge_trj(topology, trj_list, finaltrj, remove=False)
        
        return os.path.abspath(finaltrj)



def calc_avg(value_list):
    sum = 0
    l = len(value_list)

    for v in value_list:
        sum += v

    avg = sum / l

    return avg

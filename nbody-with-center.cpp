#include<iostream>
#include<fstream>
#include<unistd.h>
#include<sys/stat.h>
#include<particle_simulator.hpp>

//#include <bits/stdc++.h>
constexpr PS::F64 MY_PI = M_PI;

#ifdef ENABLE_PHANTOM_GRAPE_X86
#include <gp5util.h>
#endif
#ifdef ENABLE_GPU_CUDA
#define MULTI_WALK
#include"force_gpu_cuda.hpp"
#endif
#include "user-defined.hpp"

void mybarrier()
{
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
    PS::Comm::barrier();
#endif
}

const static PS::S32 n_tadds=10;
class ElapsedTimes{
public:
    PS::F64 t_domain_decomposition;
    PS::F64 t_exchange_particle;
    PS::F64 t_treeforce;
    PS::F64 t_treeforce2;
    PS::F64 t_misc;
    PS::F64 t_total;
    PS::F64 t_offset;
    PS::F64 t_additional[n_tadds];
    PS::S32 steps;
    void clear()
    {
	t_offset = PS::GetWtime();
	t_total=0;
	t_domain_decomposition=0;
	t_exchange_particle =0;
	t_treeforce=0;
	t_treeforce2=0;
	t_misc=0;
	steps=0;
	for(int i=0;i<n_tadds; i++)t_additional[i]=0;
    }
    void measure(int start,
		 int mode)
    {
	static PS::F64 tstart;
	static PS::S32 mode_called;
	if (start){
	    tstart = PS::GetWtime();
	    mode_called = mode;
	}else{
	    auto dt = PS::GetWtime() - tstart;
	    if (mode != mode_called){
		fprintf(stdout,
			"Invalid call to ElapsedTime#meadure:modes= %d  %d\n",
			mode_called, mode);
	    }
	    if (mode == 0){
		t_domain_decomposition += dt;
	    }else if (mode == 1){
		t_exchange_particle += dt;
	    }else if (mode == 2){
		t_treeforce+= dt;
	    }else if (mode == 3){ 
		t_treeforce2+= dt;
	    }else if (mode > 3){ 
		//		fprintf(stderr, "mode=%d dt=%e\n", mode, dt);
		t_additional[mode-4]+=dt;
	    }
	}
    }
    void inc_count()
    {
	steps++;
    }
    void print(FILE* f,
	       PS::Count_t ncep,
	       PS::Count_t ncsp)
    {
	if(PS::Comm::getRank() == 0){
	    t_total = PS::GetWtime() - t_offset;
	    t_misc = t_total - t_domain_decomposition - t_exchange_particle
		- t_treeforce - t_treeforce2;
	    if (steps==0)return;
	    fprintf(f," Total time=%.4g steps=%d, t=%.4g dd=%.4g exch=%.4g force=%.4g,%.4g misc=%.4g",
		    t_total, steps, t_total/steps,
		    t_domain_decomposition/steps,
		    t_exchange_particle/steps,
		    t_treeforce/steps,
		    t_treeforce2/steps,
		    t_misc/steps);
	    fprintf(f," ints= %.4g %.4g %.4g\n",
		    ((double)ncep)/steps,
		    ((double)ncsp)/steps,
		    ((double)(ncep+ncsp))/steps);
	    for(int i=0;i<4;i++){
		fprintf(f," ta[%d]=%.4g", i, t_additional[i]/steps);
	    }
	    fprintf(f,"\n");
	    
	}
	
    }
};

		
template<class Tpsys>
void generate_initial_cold_disk(PS::S64 nx,
				PS::S64 ny,
				int np,
				int myid,
				PS::F64 width,
				PS::F64 mass,
				Tpsys & psys,
				PS::S32 & n_loc)
{
    PS::S64 mymin = (nx+0.5)/np*myid;
    PS::S64 mymax = (nx+0.5)/np*(myid+1);
    //    fprintf(stderr, "initial disk nx=%lld ny=%lld width=%e  np=%d, id=%d, min=%lld, max==%lld\n",
    //	    nx,ny, width,  np, myid, mymin, mymax);
    if (mymax  > nx) mymax = nx;
    auto nxlocal = mymax - mymin;
    n_loc = nxlocal * ny;
    psys.setNumberOfParticleLocal(n_loc);
    PS::S64 id = mymin * ny;
    PS::S64 i=0;
    for(PS::S32 ix = 0; ix < nxlocal; ix++){
	PS::F64 theta = (ix + mymin+0.5)*2*M_PI/nx - M_PI;
	PS::F64 x0 = cos(theta);
	PS::F64 y0 = sin(theta);
	PS::F64 dr;
	if (ny > 1){
	    dr=  width /(ny -1);
	}else{
	    dr = 0;
	}
	//	fprintf(stderr, "ix  x0, y0, theta, width, dr = %d %e %e %e %e %e\n",
	//		ix, x0, y0, theta, width, dr);
	for(PS::S32 iy = 0; iy < ny; iy++){
	    PS::F64 rscale = (iy* dr) + 1.0 - width/2;
	    auto vscale = 1.0/sqrt(rscale);
	    psys[i].pos_car  = PS::F64vec(x0*rscale, y0*rscale, 0.0);
	    psys[i].vel  = PS::F64vec(-y0*vscale, x0*vscale, 0.0);
	    psys[i].mass  = mass;
	    psys[i].id   = id;
	    i++;
	    id++;
	    
	}
    }
}

	
    
    

void makeColdUniformSphere(const PS::F64 mass_glb,
                           const PS::S64 n_glb,
                           const PS::S64 n_loc,
                           PS::F64 *& mass,
                           PS::F64vec *& pos,
                           PS::F64vec *& vel,
                           const PS::F64 eng = -0.25,
                           const PS::S32 seed = 0) {
    
    assert(eng < 0.0);
    {
        PS::MTTS mt;
        mt.init_genrand(0);
        for(PS::S32 i = 0; i < n_loc; i++){
            mass[i] = mass_glb / n_glb;
            const PS::F64 radius = 3.0;
            do {
                pos[i][0] = (2. * mt.genrand_res53() - 1.) * radius;
                pos[i][1] = (2. * mt.genrand_res53() - 1.) * radius;
                pos[i][2] = (2. * mt.genrand_res53() - 1.) * radius;
            }while(pos[i] * pos[i] >= radius * radius);
            vel[i][0] = 0.0;
            vel[i][1] = 0.0;
            vel[i][2] = 0.0;
        }
    }

    PS::F64vec cm_pos  = 0.0;
    PS::F64vec cm_vel  = 0.0;
    PS::F64    cm_mass = 0.0;
    for(PS::S32 i = 0; i < n_loc; i++){
        cm_pos  += mass[i] * pos[i];
        cm_vel  += mass[i] * vel[i];
        cm_mass += mass[i];
    }
    cm_pos /= cm_mass;
    cm_vel /= cm_mass;
    for(PS::S32 i = 0; i < n_loc; i++){
        pos[i] -= cm_pos;
        vel[i] -= cm_vel;
    }
}

template<class Tpsys>
void setParticlesColdUniformSphere(Tpsys & psys,
                                   const PS::S32 n_glb,
                                   PS::S32 & n_loc) {

    n_loc = n_glb; 
    psys.setNumberOfParticleLocal(n_loc);

    PS::F64    * mass = new PS::F64[n_loc];
    PS::F64vec * pos  = new PS::F64vec[n_loc];
    PS::F64vec * vel  = new PS::F64vec[n_loc];
    const PS::F64 m_tot = 1.0;
    const PS::F64 eng   = -0.25;
    makeColdUniformSphere(m_tot, n_glb, n_loc, mass, pos, vel, eng);
    for(PS::S32 i = 0; i < n_loc; i++){
        psys[i].mass = mass[i];
        psys[i].pos  = pos[i];
        psys[i].vel  = vel[i];
        psys[i].id   = i;
    }
    delete [] mass;
    delete [] pos;
    delete [] vel;
}

template<class Tpsys>
void kick(Tpsys & system,
          const PS::F64 dt,
	  PS::F64 & elost_local) {
    PS::S32 n = system.getNumberOfParticleLocal();
    //    std::cerr << "Kick\n";
#pragma omp parallel for reduction(+:elost_local)
    for(PS::S32 i = 0; i < n; i++) {
	// 	std::cerr << i << " "<< system[i].acc <<
	//	    " "<< system[i].pot << " " <<system[i].pos_car << "\n" ;
        system[i].vel  += system[i].acc * dt;
	elost_local += dt* system[i].dedtdumper* system[i].mass;
	//	fprintf(stderr, "dedtdumper[%d]=%e elost=%g\n", i, system[i].dedtdumper, elost_local);
    }
}

template<class Tpsys>
bool remove_particles(Tpsys & system,
		      PS::F64 tsys,
		      PS::F64 r_inner_cutoff,
		      PS::F64 & elost_local,
		      std::ofstream & fout_remove)
{
    PS::F64 rcut2 = r_inner_cutoff*r_inner_cutoff;
    PS::S32 n = system.getNumberOfParticleLocal();
    PS::S32 removeindex[n];
    PS::S32 nremove=0;
    PS::F64 deremove = 0;
    for(PS::S32 i = 0; i < n; i++) {
	if (system[i].pos_car*system[i].pos_car < rcut2){
	    removeindex[nremove]=i;
	    nremove ++;
	    deremove += system[i].mass*(system[i].vel * system[i].vel*0.5
				       + system[i].pot
				       + system[i].pot_center);
	    fout_remove << "Remove t "<<tsys << " id= " << system[i].id <<" "
			<< system[i].mass << " "
			<< system[i].pos_car << " "
			<< system[i].vel << "\n";
	    
		
	}
    }
    if (nremove > 0){
	system.removeParticle(removeindex, nremove);
	elost_local += deremove;
    }
    return (nremove > 0);
    
}

template<class Tpsys>
void dtoc(Tpsys & system)
{
    PS::S32 n = system.getNumberOfParticleLocal();
#pragma omp parallel for
    for(PS::S32 i = 0; i < n; i++) {
        system[i].dtoc();
    }
}

template<class Tpsys>
void add_central_gravity(Tpsys & system)
{
    PS::S32 n = system.getNumberOfParticleLocal();
#pragma omp parallel for
    for(PS::S32 i = 0; i < n; i++) {
	system[i].add_central_gravity();
    }
}

template<class Tpsys>
void dump_system(Tpsys & system)
{
    PS::S32 n = system.getNumberOfParticleLocal();
    fprintf(stderr,"Dumping all particles\n");
    for(PS::S32 i = 0; i < n; i++) {
	system[i].dumpparticle();
    }
}

template<class Tpsys>
void drift(Tpsys & system,
           const PS::F64 dt) {
    PS::S32 n = system.getNumberOfParticleLocal();
#pragma omp parallel for
    for(PS::S32 i = 0; i < n; i++) {
        system[i].pos_car  += system[i].vel * dt;
    }
}

template<class Tpsys>
void calcEnergy(const Tpsys & system,
                PS::F64 & etot,
                PS::F64 & ekin,
                PS::F64 & epot,
		PS::F64 & elost,
		PS::F64 elost_local,
                const bool clear=true){
    if(clear){
        etot = ekin = epot =  elost = 0.0;
    }
    PS::F64 etot_loc = 0.0;
    PS::F64 ekin_loc = 0.0;
    PS::F64 epot_loc = 0.0;
    PS::F64 epot_center_loc = 0.0;
    const PS::S32 nbody = system.getNumberOfParticleLocal();
    for(PS::S32 i = 0; i < nbody; i++){
        ekin_loc += system[i].mass * system[i].vel * system[i].vel;
        epot_loc += system[i].mass * system[i].pot;
        epot_center_loc += system[i].mass * system[i].pot_center;
    }
    ekin_loc *= 0.5;
    epot_loc *= 0.5;
    epot_loc += epot_center_loc;
    etot_loc  = ekin_loc + epot_loc + elost_local;
    etot = PS::Comm::getSum(etot_loc);
    epot = PS::Comm::getSum(epot_loc);
    ekin = PS::Comm::getSum(ekin_loc);
    elost = PS::Comm::getSum(elost_local);
    //    fprintf(stderr, "ekin, epot, etot= %g %g %g\n", ekin, epot, etot);
}

void printHelp() {
    std::cerr<<"b: show time profile"<<std::endl;
    std::cerr<<"D: dt_snap (default: 1.0)"<<std::endl;
    std::cerr<<"d: dt_diag (default: 1.0 / 8.0)"<<std::endl;
    std::cerr<<"E: set exchange let mode (0,1,2, default=0)"<<std::endl;
    std::cerr<<"e: eta (default: 0.0, set if K option>0)"<<std::endl;
    std::cerr<<"g: inner radius to remove particles (default: 0.0)"<<std::endl;
    std::cerr<<"I: Initial disk model (nx, ny, width, mass. default=1000,20,0.1,5r-8)"<<std::endl;
    std::cerr<<"i: dir name of input (default:none)"<<std::endl;
    std::cerr<<"K: ndt for bound (default: 32)"<<std::endl;
    std::cerr<<"k: kappa (default: 0.0, set if K option>0)"<<std::endl;
    std::cerr<<"l: n_leaf_limit (default: 8)"<<std::endl;
    std::cerr<<"M: central mass (default: 1)"<<std::endl;
    std::cerr<<"N: n_tot (default: 1024)"<<std::endl;
    std::cerr<<"n: n_group_limit (default: 64)"<<std::endl;
    std::cerr<<"o: dir name of output (default: ./result)"<<std::endl;
    std::cerr<<"p: eps (default: 0)"<<std::endl;
    std::cerr<<"R: number of divisions in radial direction   (default: 1)"<<std::endl;
    std::cerr<<"r: rcoll radius of particles (default: 0.0)"<<std::endl;
    std::cerr<<"S: replulsion coef (default: 0.5)"<<std::endl;
    std::cerr<<"s: time_step (default: 1.0 / 512.0)"<<std::endl;
    std::cerr<<"T: time_end (default: 10.0)"<<std::endl;
    std::cerr<<"t: theta (default: 0.5)"<<std::endl;
    std::cerr<<"h: help"<<std::endl;
}

void makeOutputDirectory(char * dir_name) {
    struct stat st;
    PS::S32 ret;
    if (dir_name[0]==0) return;
    if (PS::Comm::getRank() == 0) {
        if (stat(dir_name, &st) != 0) {
            ret = mkdir(dir_name, 0777);
        } else {
            ret = 0; // the directory named dir_name already exists.
        }
    } 
    PS::Comm::broadcast(&ret, 1);
    if (ret == 0) {
        if (PS::Comm::getRank() == 0)
            fprintf(stderr, "Directory \"%s\" is successfully made.\n", dir_name);
    } else {
        if (PS::Comm::getRank() == 0)
            fprintf(stderr, "Directory %s fails to be made.\n", dir_name);
        PS::Abort();
    }
}


void set_coeffs( PS::F64 *kappa,
		 PS::F64 *eta,
		 PS::F32 repl_coef,
		 PS::F32 period)
{
    *eta = -4*log(repl_coef)/period;
    *kappa = 4*M_PI*M_PI/(period*period) + (*eta)*(*eta)/4.0;
}

template<class Tpsys>
void printtimeprofile(PS::TimeProfile tp,
		      PS::TimeProfile tps,
		      PS::TimeProfile tpd,
		      FILE* f,
		      Tpsys &system, // particle system
		      int steps)
{
    if (steps==0)steps=1;
    auto nsend = system.getNumberOfParticleSendGlobal();
    PS::F64  nslocal = (PS::F64)system.getNumberOfParticleSendLocal();
    auto nsmax = PS::Comm::getMaxValue(nslocal);
    auto epepmax = PS::Comm::getMaxValue(tps.exchange_particle__exchange_particle);
    auto epepmin = PS::Comm::getMinValue(tps.exchange_particle__exchange_particle);
    auto epfmax = PS::Comm::getMaxValue(tps.exchange_particle__find_particle);
    auto epfmin = PS::Comm::getMinValue(tps.exchange_particle__find_particle);
    auto epepmax1=PS::Comm::getMaxValue(tps.exchange_particle__exchange_particle_1);
    auto epepmin1=PS::Comm::getMinValue(tps.exchange_particle__exchange_particle_1);
    auto epepmax2=PS::Comm::getMaxValue(tps.exchange_particle__exchange_particle_2);
    auto epepmax3=PS::Comm::getMaxValue(tps.exchange_particle__exchange_particle_3);
    auto el1max=PS::Comm::getMaxValue(tp.exchange_LET_1st);
    auto el1min=PS::Comm::getMinValue(tp.exchange_LET_1st);
    auto el1topmax=PS::Comm::getMaxValue(tp.make_LET_1st__exchange_top_moment);
    auto el1tfp=PS::Comm::getMaxValue(tp.make_LET_1st__find_particle);
    auto el1tnmax=PS::Comm::getMaxValue(tp.make_LET_1st__exchange_n);
    auto el1ticp=PS::Comm::getMaxValue(tp.exchange_LET_1st__icomm_ptcl);
    auto cfmax=PS::Comm::getMaxValue(tp.calc_force);
    auto cfwt=PS::Comm::getMaxValue(tp.calc_force__core__walk_tree);
    auto cfmkipg=PS::Comm::getMaxValue(tp.calc_force__make_ipgroup);
    auto cfc=PS::Comm::getMaxValue(tp.calc_force__core);
    auto cs=PS::Comm::getMaxValue(tpd.collect_sample_particle);
    auto dd=PS::Comm::getMaxValue(tpd.decompose_domain);
    auto ddgp=PS::Comm::getMaxValue(tpd.decompose_domain__gather_particle);
    auto dds1=PS::Comm::getMaxValue(tpd.decompose_domain__sort_particle_1st);
    auto dds2=PS::Comm::getMaxValue(tpd.decompose_domain__sort_particle_2nd);
    auto dds3=PS::Comm::getMaxValue(tpd.decompose_domain__sort_particle_3rd);
    auto dded=PS::Comm::getMaxValue(tpd.decompose_domain__exchange_pos_domain);
    auto ddb=PS::Comm::getMaxValue(tpd.decompose_domain__barrier);
    auto ddp=PS::Comm::getMaxValue(tpd.decompose_domain__postprocess);
    auto ep=PS::Comm::getMaxValue(tps.exchange_particle);
    auto ml=PS::Comm::getMaxValue(tp.make_local_tree);
    auto mg=PS::Comm::getMaxValue(tp.make_global_tree);
    auto cf=PS::Comm::getMaxValue(tp.calc_force);
    auto cml=PS::Comm::getMaxValue(tp.calc_moment_local_tree);
    auto cmg=PS::Comm::getMaxValue(tp.calc_moment_global_tree);
    auto ml1=PS::Comm::getMaxValue(tp.make_LET_1st);
    auto ml2=PS::Comm::getMaxValue(tp.make_LET_2nd);
    auto el1=PS::Comm::getMaxValue(tp.exchange_LET_1st);
    auto el2=PS::Comm::getMaxValue(tp.exchange_LET_2nd);
    
    if(PS::Comm::getRank() == 0){
	fprintf(f, " cs: %.4g dd: %.4g ep:%.4g ml:%.4g mg:%.4g cf:%.4g cml:%.4g cmg:%.4g lets: %.4g %.4g %.4g %.4g\n",
		cs/steps,
		dd/steps,
		ep/steps,
		ml/steps,
		mg/steps,
		cf/steps,
		cml/steps,
		cmg/steps,
		ml1/steps,
		ml2/steps,
		el1/steps,
		el2/steps);
	auto alltp = tp+tps+tpd;
	fprintf(f, " epf=%.4g %.4g %.4g epe=%.4g %.4g %.4g  mkl=%.4g msl=%.4g lcl=%.4g msg=%.4g lcg=%.4g ",
		alltp.exchange_particle__find_particle/steps,
		epfmax/steps,
		epfmin/steps,
		alltp.exchange_particle__exchange_particle/steps,
		epepmax/steps,
		epepmin/steps,
		tp.morton_key_local_tree/steps,
		alltp.morton_sort_local_tree/steps,
		alltp.link_cell_local_tree/steps,
		alltp.morton_sort_global_tree/steps,
		alltp.link_cell_global_tree/steps);
	fprintf(f, "ns=%.4g nsmax=%.4g\n",
		(PS::F64) nsend/steps,
		(PS::F64) nsmax/steps);
#if 1	 //details for domain decomposition
	fprintf(f, " ddgp=%.4g s1=%.4g s2=%.4g s3=%.4g b=%.4g ed=%.4g p=%.4g",
		ddgp/steps,
		dds1/steps,
		dds2/steps,
		dds3/steps,
		ddb/steps,
		dded/steps,
		ddp/steps);
#endif	
#if 1	 //details for exchange particles
	fprintf(f, " epeptimes=%.4g %.4g %.4g %.4g",
		epepmax1/steps,
		epepmin1/steps,
		epepmax2/steps,
		epepmax3/steps);
#endif	
	fprintf(f, " el1=%.4g %.4g tm=%.4g fp=%.4g en=%.4g commp=%.4g",
		el1max/steps,
		el1min/steps,
		el1topmax/steps,
		el1tfp/steps,
		el1tnmax/steps,
		el1ticp/steps);

	fprintf(f,"\n");
	fprintf(f, " cfmax=%.4g wt=%.4g ipg=%.4g core=%.4g",
		cfmax/steps,
		cfwt/steps,
		cfmkipg/steps,
		cfc/steps	);
	fprintf(f,"\n");

	
		
    }
    
}
    
PS::F64  FPGrav::eps = 0;
PS::F64    FPGrav::rcoll = 0;
PS::F64    FPGrav::kappa = 0;
PS::F64    FPGrav::eta = 0;
PS::F64    FPGrav::mass_center = 1.0;


#ifdef QUAD
using SPJ_t    = MySPJQuadrupole;
using Moment_t = MyMomentQuadrupole;
using CalcForceSp = CalcForceSpQuad<FPGrav, SPJ_t, ForceGrav>;
#else
using SPJ_t    = MySPJMonopole;
using Moment_t = MyMomentMonopole;
using CalcForceSp = CalcForceSpMono<FPGrav, SPJ_t, ForceGrav>;
#endif
    
using MY_SEARCH_MODE = PS::SEARCH_MODE_LONG_SCATTER;
//using MY_SEARCH_MODE = PS::SEARCH_MODE_LONG_SYMMETRY;
// the use of SEARCH_MODE_LONG_SYMMETRY for test purpose only....
//using Tree_t = PS::TreeForForce<MY_SEARCH_MODE, ForceGrav, FPGrav, EPJGrav, Moment_t, Moment_t, SPJ_t, PS::CALC_DISTANCE_TYPE_NORMAL>;
using Tree_t = PS::TreeForForce<MY_SEARCH_MODE, ForceGrav, FPGrav, EPJGrav, Moment_t, Moment_t, SPJ_t, PS::CALC_DISTANCE_TYPE_NEAREST_X>;

int main(int argc, char *argv[]) {
    std::cout<<std::setprecision(15);
    std::cerr<<std::setprecision(15);

    
    PS::Initialize(argc, argv);
    PS::F32 theta = 0.5;
    PS::S32 n_leaf_limit = 8;
    PS::S32 n_group_limit = 64;
    PS::F32 time_end = 10.0;
    PS::F32 dt = 1.0 / 1024.0;
    PS::F32 dt_diag = 1.0 / 8.0;
    PS::F32 dt_snap = 1.0;
    PS::F32 ndtbound = 32;
    PS::F32 repl_coef = 0.5;
    PS::S32 ny_for_dd = 1;
    PS::S32 show_time_profile = 0;
    char dir_name[1024];
    PS::S32 makeoutput=0;
    char in_name[1024];
    PS::S64 n_tot = 1024;
    PS::S32 c;
    PS::S32 exchange_let_mode = 0;
    PS::S64 nxinit = 1000;
    PS::S64 nyinit = 20;
    PS::F64 widthinit = 0.1;
    PS::F64 massinit = 5e-8;
    PS::F64 r_inner_cutoff = 0.0;
    //    strncpy(dir_name,"./result", 1000);
    dir_name[0]=0;
    in_name[0]=0;
    opterr = 0;
    while((c=getopt(argc,argv,"i:r:k:e:g:p:o:d:D:t:T:l:M:n:hs:K:S:R:bE:I:")) != -1){
        switch(c){
	    case 'i':
		strncpy(in_name,optarg,1000);
		break;
	    case 'o':
		//            sprintf(dir_name,optarg);
		strncpy(dir_name,optarg,1000);
		makeoutput = 1;
		break;
	    case 'M':
		FPGrav::mass_center	= atof(optarg);
		break;
	    // case 'k':
	    // 	FPGrav::kappa = atof(optarg);
	    // 	break;
	    case 'p':
		FPGrav::eps  = atof(optarg);
		break;
	    case 'e':
		FPGrav::eta  = atof(optarg);
		break;
	    case 'g':
		r_inner_cutoff  = atof(optarg);
		break;
	    case 'r':
		FPGrav::rcoll  = atof(optarg);
		break;
	    case 't':
		theta = atof(optarg);
		std::cerr << "theta =" << theta << std::endl;
		break;
	    case 'T':
		time_end = atof(optarg);
		break;
	    case 's':
		dt = atof(optarg);
		break;
	    case 'd':
		dt_diag = atof(optarg);
		break;
	    case 'D':
		dt_snap = atof(optarg);
		break;
	    case 'l':
		n_leaf_limit = atoi(optarg);
		break;
	    case 'n':
		n_group_limit = atoi(optarg);
		break;
	    case 'K':
		ndtbound = atof(optarg);
		break;
	    case 'R':
		ny_for_dd = atoi(optarg);
		break;
	    case 'S':
		repl_coef = atof(optarg);
		break;
	    case 'b':
		show_time_profile = 1;
		break;
	    case 'E':
		exchange_let_mode = atoi(optarg);
		break;
	    case 'I':
		sscanf(optarg,"%lld,%lld,%lf,%lf", &nxinit, &nyinit, &widthinit, &massinit);
		break;
	    case 'h':
		if(PS::Comm::getRank() == 0) {
		    printHelp();
		}
		PS::Finalize();
		return 0;
	    default:
		if(PS::Comm::getRank() == 0) {
		    std::cerr<<"No such option! Available options are here."<<std::endl;
		    printHelp();
		}
		PS::Abort();
        }
    }
    if (ndtbound > 0.0 )set_coeffs( &(FPGrav::kappa), &(FPGrav::eta),
				    repl_coef, ndtbound * dt);
    
    if(PS::Comm::getRank() == 0) {
	std::cerr << "central_mass =" << FPGrav::mass_center<< std::endl;
	std::cerr << "kappa =" << FPGrav::kappa<< std::endl;
	std::cerr << "theta =" << theta << std::endl;
	std::cerr << "eps =" << FPGrav::eps << std::endl;
	std::cerr << "eta =" << FPGrav::eta << std::endl;
	std::cerr << "repl =" << repl_coef << std::endl;
	std::cerr << "Tcoef =" << ndtbound<< std::endl;
	std::cerr << "rcoll =" << FPGrav::rcoll << std::endl;
	std::cerr << "time_end = " << time_end << std::endl;
	std::cerr << "time_step = " << dt << std::endl;
	std::cerr << "dt_diag = " << dt_diag << std::endl;
	std::cerr << "dt_snap = " << dt_snap << std::endl;
	std::cerr << "n_leaf_limit = " << n_leaf_limit << std::endl;
	std::cerr << "n_group_limit = " << n_group_limit << std::endl;
	//	std::cerr << "n_tot = " << n_tot << std::endl;
	std::cerr << "ny_for_dd = " << ny_for_dd << std::endl;
	std::cerr << "exchange_let_mode = " << exchange_let_mode<< std::endl;
	std::cerr << "nxinit = " << nxinit<< std::endl;
	std::cerr << "nyinit = " << nyinit<< std::endl;
	std::cerr << "widthinit = " << widthinit<< std::endl;
	std::cerr << "massinit = " << massinit<< std::endl;
	std::cerr << "r_inner_cutoff = " << r_inner_cutoff<< std::endl;
    }
    makeOutputDirectory(dir_name);

    std::ofstream fout_eng;
    std::ofstream fout_remove;

    if(PS::Comm::getRank() == 0) {
        char sout_de[1024];
	sprintf(sout_de, "%s/t-de.dat", dir_name);
        fout_eng.open(sout_de);
        fprintf(stdout, "This is a sample program of N-body simulation on FDPS!\n");
        fprintf(stdout, "Number of processes: %d\n", PS::Comm::getNumberOfProc());
        fprintf(stdout, "Number of threads per process: %d\n", PS::Comm::getNumberOfThread());
    }

    {
        char sout_remove[1024];
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
	sprintf(sout_remove, "%s/remove%-4d.dat", dir_name,PS::Comm::getRank() );
#else
	sprintf(sout_remove, "%s/remove.dat", dir_name);
#endif	
        fout_remove.open(sout_remove);
    }

    PS::ParticleSystem<FPGrav> system_grav;
    system_grav.initialize();
    system_grav.clearTimeProfile();
    system_grav.clearCounterAll();
    PS::S32 n_loc    = 0;
    PS::F32 time_sys = 0.0;
    //    fprintf(stderr,"process %d enter initial condition\n",  PS::Comm::getRank());
    if (in_name[0]==0){
	generate_initial_cold_disk(nxinit, nyinit,
				   PS::Comm::getNumberOfProc(),
				   PS::Comm::getRank(),
				   widthinit,
				   massinit,
				   system_grav, n_loc);
	
    }else{
	FileHeader header;
	//	fprintf(stderr,"reading in_name = %s\n",in_name);
	system_grav.readParticleAscii(in_name, header);
	n_loc = system_grav.getNumberOfParticleLocal();
	n_tot = PS::Comm::getSum(n_loc);
    }
    
    //    fprintf(stderr,"process%d np=%d nt=%d\n",  PS::Comm::getRank(), n_loc, n_tot);
    const PS::F32 coef_ema = 0.3;
    PS::DomainInfo dinfo;
    dinfo.initialize(coef_ema);
    dinfo.clearTimeProfile();
    if (PS::Comm::getNumberOfProc() > 1){
	int nproc = PS::Comm::getNumberOfProc();
	int ny = ny_for_dd;
	int nx = nproc/ny;
	while (nx*ny != nproc){
	    nx++;
	    ny = nproc/nx;
	}
	if (PS::Comm::getRank()==0){
	    fprintf(stderr, "np = %d, nx=%d, ny=%d\n", nproc, nx, ny);
	}
	dinfo.setDomain(nx, ny);
    }
    dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_OPEN);
    //    std::cerr << "MY_PI = " << MY_PI << "\n";
    dinfo.setPosRootDomainX(-MY_PI, MY_PI);
    
    dtoc(system_grav);
    PS::F64 time0= PS::GetWtime();
    dinfo.decomposeDomainAll(system_grav);

    n_loc = system_grav.getNumberOfParticleLocal();
    for(auto i=0; i<n_loc; i++){
	if(!dinfo.getPosRootDomain().contained(system_grav[i].getPos())){
	    std::cerr<<"n_loc= "<<n_loc<<std::endl;
	    std::cerr<<"dinfo.getPosRootDomain()= "<<dinfo.getPosRootDomain()<<std::endl;
	    std::cerr<<"i= "<<i<<" system_grav[i].getPos()= "<<system_grav[i].getPos()<<std::endl;
	    std::cerr<<"i= "<<i<<" system_grav[i].pos_car= "<<system_grav[i].pos_car<<std::endl;
	    std::cerr<<"i= "<<i<<" system_grav[i].pos= "<<system_grav[i].pos<<std::endl;
	}
	assert(dinfo.getPosRootDomain().contained(system_grav[i].getPos()));
    }
    
    system_grav.exchangeParticle(dinfo);
    n_loc = system_grav.getNumberOfParticleLocal();

    
#ifdef ENABLE_PHANTOM_GRAPE_X86
    g5_open();
    g5_set_eps_to_all(FPGrav::eps);
#endif
    


    Tree_t tree_grav;
    
    tree_grav.initialize(n_tot, theta, n_leaf_limit, n_group_limit);
    auto exchange_text = "A2A";
    if (exchange_let_mode == 0){
	tree_grav.setExchangeLETMode(PS::EXCHANGE_LET_A2A);
    }else if (exchange_let_mode == 1){
	tree_grav.setExchangeLETMode(PS::EXCHANGE_LET_P2P_EXACT);
	exchange_text = "P2P_EXACT";
    }else if (exchange_let_mode == 2){
	tree_grav.setExchangeLETMode(PS::EXCHANGE_LET_P2P_FAST);
	exchange_text = "P2P_FAST";
    }
    if(PS::Comm::getRank() == 0){
	fprintf(stdout, " EXCHANGE_LET_MODE = %s\n", exchange_text);
    }
    tree_grav.clearTimeProfile();
    tree_grav.clearNumberOfInteraction();
    
#ifdef MULTI_WALK
    const PS::S32 n_walk_limit = 200;
    const PS::S32 tag_max = 1;
    tree_grav.calcForceAllAndWriteBackMultiWalk(DispatchKernelWithSP,
                                                RetrieveKernel,
                                                tag_max,
                                                system_grav,
                                                dinfo,
                                                n_walk_limit);
#else
    tree_grav.calcForceAllAndWriteBack(CalcForceEp<EPJGrav>,
                                       CalcForceSp(),
                                       system_grav,
                                       dinfo);
    //				       true, PS::MAKE_LIST_FOR_REUSE );
#endif
    add_central_gravity(system_grav);
    //    dump_system(system_grav);
    PS::F64 Epot0, Ekin0, Etot0, Epot1, Ekin1, Etot1, Elost0, Elost1,Elost_local;
    calcEnergy(system_grav, Etot0, Ekin0, Epot0, Elost0, 0.0);
    Elost_local=0.0;
    PS::F64 time_diag = 0.0;
    PS::F64 time_snap = 0.0;
    PS::S64 n_loop = 0;
    PS::S64 n_loop_diag = 0;
    PS::S32 id_snap = 0;
    ElapsedTimes et;
    et.clear();
    while(time_sys < time_end){
	et.measure(1,6);
	
        if( (time_sys >= time_snap) || ( (time_sys + dt) - time_snap ) > (time_snap - time_sys) ){
            char filename[256];
            sprintf(filename, "%s/%04d.dat", dir_name, id_snap++);
            FileHeader header;
            header.time   = time_sys;
            header.n_body = system_grav.getNumberOfParticleGlobal();
            if (makeoutput) system_grav.writeParticleAscii(filename, header);
            time_snap += dt_snap;
        }

        calcEnergy(system_grav, Etot1, Ekin1, Epot1, Elost1, Elost_local);
	if (show_time_profile)   mybarrier();
        
	if( (time_sys >= time_diag) || ( (time_sys + dt) - time_diag ) > (time_diag - time_sys) ){
	    PS::Count_t ncepglobal = tree_grav.getNumberOfInteractionEPEPGlobal();
	    PS::Count_t ncspglobal = tree_grav.getNumberOfInteractionEPSPGlobal();
	    PS::TimeProfile tp = tree_grav.getTimeProfile();
	    PS::TimeProfile tps = system_grav.getTimeProfile();
	    PS::TimeProfile tpd = dinfo.getTimeProfile();
	    
	    if(PS::Comm::getRank() == 0){
                fout_eng << time_sys << "   " << (Etot1 - Etot0) / Etot0 <<
		    " " << Etot1 << " " << Ekin1 << " " << Epot1 <<std::endl;
                fprintf(stdout, "time: %10.7f etot,k,p,l=%e %e %e %e energy error: %+e\n",
                        time_sys, Etot1, Ekin1, Epot1, Elost1,
			(Etot1 - Etot0) / Etot0);

            }
	    if (show_time_profile){
		printtimeprofile(tp, tps, tpd, stdout, system_grav, et.steps);
		et.print(stdout,ncepglobal, ncspglobal  );
		//  fprintf(stdout,"Wall clock time=%10.5f\n", PS::GetWtime()-time0);
	    }
	    if(PS::Comm::getRank() == 0 && n_loop_diag %32== 0) fflush(stdout);
	    n_loop_diag ++;
	    tree_grav.clearTimeProfile();
	    dinfo.clearTimeProfile();
	    system_grav.clearTimeProfile();
	    system_grav.clearCounterAll();
	    tree_grav.clearNumberOfInteraction();
	    et.clear();
	    time_diag += dt_diag;
	    
        }
        
	et.measure(0,6);
	et.measure(1,5);
        
        kick(system_grav, dt * 0.5, Elost_local);
        
        time_sys += dt;
        drift(system_grav, dt);
        dtoc(system_grav);
	if (show_time_profile)   mybarrier();
	et.measure(0,5);
        
	et.measure(1,0);
        if(n_loop % 4 == 0){
            dinfo.decomposeDomainAll(system_grav);
        }
        
	et.measure(0,0);
	et.measure(1,7);
	if (show_time_profile)   mybarrier();
	et.measure(0,7);
	et.measure(1,1);
        system_grav.exchangeParticle(dinfo);
	et.measure(0,1);
	
	et.measure(1,2);
#ifdef MULTI_WALK
        tree_grav.calcForceAllAndWriteBackMultiWalk(DispatchKernelWithSP,
                                                    RetrieveKernel,
                                                    tag_max,
                                                    system_grav,
                                                    dinfo,
                                                    n_walk_limit,
                                                    true);
#else
        tree_grav.calcForceAllAndWriteBack(CalcForceEp<EPJGrav>,
                                           CalcForceSp(),
                                           system_grav,
                                           dinfo);
	//		   true, PS::MAKE_LIST_FOR_REUSE );
	
#endif
	et.measure(0,2);
	et.measure(1,3);
	if (show_time_profile) mybarrier();
	et.measure(0,3);
	et.measure(1,4);
	add_central_gravity(system_grav);
        kick(system_grav, dt * 0.5, Elost_local);
	//	dump_system(system_grav);
	remove_particles(system_grav, time_sys,  r_inner_cutoff, Elost_local,
			 fout_remove);
	
        
	if (show_time_profile) mybarrier();
	et.measure(0,4);
        n_loop++;
	et.inc_count();
    }
    
#ifdef ENABLE_PHANTOM_GRAPE_X86
    g5_close();
#endif

    PS::Finalize();
    return 0;
}


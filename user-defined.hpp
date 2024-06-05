#pragma once
//constexpr PS::F64 MY_PI = M_PI;
constexpr PS::F64 SMALL = 1e-15;
class FileHeader{
public:
    PS::S64 n_body;
    PS::F64 time;
    PS::S32 readAscii(FILE * fp) {
	int i= fscanf(fp, "%lf\n", &time);
        i+= fscanf(fp, "%lld\n", &n_body);
        return n_body;
    }
    void writeAscii(FILE* fp) const {
        fprintf(fp, "%e\n", time);
        fprintf(fp, "%lld\n", n_body);
    }
};

class ForceGrav{
    
    
public:
    PS::F64vec acc;
    PS::F64    pot;    
    PS::F64    dedtdumper;
    
    
    void clear() {
        acc = 0.0;
        pot = 0.0;
	dedtdumper = 0.0;
    }
};

    
class FPGrav{
    
public:
    PS::S64    id;
    PS::F64    mass;
    PS::F64vec pos;
    PS::F64vec pos_car;
    PS::F64vec vel;
    PS::F64vec acc;
    PS::F64    dedtdumper;
    PS::F64    pot;    
    PS::F64    pot_center;
    static PS::F64 rcoll;
    static PS::F64 kappa;
    static PS::F64 eta;
    static PS::F64 eps;
    static PS::F64 mass_center;

    PS::F64vec getPos() const {
        return pos;
    }

    PS::F64vec getPosCar() const {
        return pos_car;
    }

    void dtoc()
    {
//	const auto pos_r   = sqrt(pos_car.x*pos_car.x + pos_car.y*pos_car.y);
	const auto pos_r   = log(pos_car.x*pos_car.x + pos_car.y*pos_car.y)*0.5+1;
	auto pos_phi = atan2(pos_car.y, pos_car.x);
	if (pos_phi >= M_PI) pos_phi -= SMALL;
	if (pos_phi <= -M_PI) pos_phi += SMALL;
	pos =  PS::F64vec(pos_phi, pos_r, pos_car.z);
    }
    // void ctod()
    // {
    // 	const auto cth = cos(pos.x);
    // 	const auto sth = sin(pos.x);
    // 	const auto r = pos.y;
    // 	const auto pos_x = r*cth;
    // 	const auto pos_y = r*sth;
    // 	pos_car= PS::F64vec(pos_x, pos_y, pos.z);
    // }
    
	
    PS::F64 getRSearch() const {
        return rcoll;
    }
    PS::F64 getCharge() const {
        return mass;
    }

    void copyFromFP(const FPGrav & fp){ 
        id = fp.id;
        mass = fp.mass;
        pos  = fp.pos;
        pos_car  = fp.pos_car;
        rcoll  = fp.rcoll;
        vel  = fp.vel;
    }

    void copyFromForce(const ForceGrav & force) {
        acc = force.acc;
        dedtdumper= force.dedtdumper;
        pot = force.pot;
    }

    void clear() {
        acc = 0.0;
        pot = 0.0;
	dedtdumper=0.0;
    }

    void dumpparticle() const {
        fprintf(stderr, "%lld\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n", 
                this->id, this->mass,
                this->pos.x, this->pos.y, this->pos.z,
                this->vel.x, this->vel.y, this->vel.z,
                this->acc.x, this->acc.y, this->acc.z,
                this->pot,   this->pot_center);
    }

    void writeAscii(FILE* fp) const {
        fprintf(fp, "%lld\t%.10g\t%.10g\t%.10g\t%.10g\t%.10g\t%.10g\t%.10g\n", 
                this->id, this->mass,
                this->pos_car.x, this->pos_car.y, this->pos_car.z,
                this->vel.x, this->vel.y, this->vel.z);
    }

    void readAscii(FILE* fp) {
        if (fscanf(fp, "%lld\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", 
		   &this->id, &this->mass,
		   &this->pos_car.x, &this->pos_car.y, &this->pos_car.z,
		   &this->vel.x, &this->vel.y, &this->vel.z) != 8){
	    fprintf(stderr,"ReadAscii failed\n");
	}
    }

    void readBinary(FILE* fp) {       }

    void add_central_gravity()
    {
	PS::F64 r2=pos_car*pos_car ;
	acc -= FPGrav::mass_center* pos_car/(r2*sqrt(r2));
	pot_center = -FPGrav::mass_center/(sqrt(r2));
	//	fprintf(stderr,"pos = %g %g %g r2= %g pot=%g pot_center= %g\n",
	//	pos.x, pos.y, pos.z, r2, pot, pot_center);
    }
    
};

class EPJGrav{
    
public:
    PS::S64    id;
    PS::F64    mass;
    PS::F64vec pos;
    PS::F64vec pos_car;
    PS::F64vec vel;
    PS::F64vec getPos() const {
        return pos;
    }

    PS::F64vec getPosCar() const {
        return pos_car;
    }

	
    PS::F64 getRSearch() const {
        return FPGrav::rcoll;
    }
    PS::F64 getCharge() const {
        return mass;
    }

    void copyFromFP(const FPGrav & fp){ 
        id = fp.id;
        mass = fp.mass;
        pos  = fp.pos;
        pos_car  = fp.pos_car;
        vel  = fp.vel;
    }

};


class MyMomentMonopole{
public:
    PS::F64 mass;
    PS::F64vec pos;
    PS::F64vec pos_car;
    PS::F64ort boundary;
    MyMomentMonopole() : mass(0.0), pos(PS::F64vec(0.0)), pos_car(PS::F64vec(0.0)), boundary(PS::F64ort(PS::F64vec(-100.0), PS::F64vec(100.0))){}
    MyMomentMonopole(const PS::F64 m, const PS::F64vec & p, const PS::F64vec & p_car, const PS::F64ort & b) : mass(m), pos(p), pos_car(p_car), boundary(b){}
    void init(){
        mass = 0.0;
        pos = 0.0;
        pos_car = 0.0;
        boundary.init();
    }
    PS::F64vec getPos() const {
        return pos;
    }
    PS::F64 getCharge() const {
	return mass;
    }
    template<class Tepj>
    void accumulateAtLeaf(const Tepj & epj){
	//	std::cerr << "leaf mass, poss= "<<epj.getCharge() << " "<< epj.getPosCar()
	//		  << " " <<epj.getPos()<<"\n";
        mass += epj.getCharge();
        pos += epj.getCharge() * epj.getPos();
        pos_car += epj.getCharge() * epj.getPosCar();
        boundary.merge(epj.getPos());
    }
    template<class Tepj>
    void accumulateAtLeaf2(const Tepj & epj){}
    void set(){
        pos = pos / mass;
        pos_car = pos_car / mass;
	//	std::cerr << "up node mass, poss= "<<mass << " "<< pos_car << " " <<pos<<"\n";
    }
    void accumulate(const MyMomentMonopole & mom){
	//	std::cerr << "node mass, poss= "<<mom.mass << " "<< mom.pos_car
	//			  << " " <<mom.pos <<"\n";
        mass += mom.mass;
        pos += mom.mass * mom.pos;
        pos_car += mom.mass * mom.pos_car;
        boundary.merge(mom.boundary);
    }
    void accumulate2(const MyMomentMonopole & mom){}
    // for DEBUG 
    void dump(std::ostream & fout = std::cout) const {
        fout<<"mass="<<mass<<std::endl;
        fout<<"pos="<<pos<<std::endl;
        fout<<"pos_car="<<pos_car<<std::endl;
    }
};

class MySPJMonopole{
public:
    PS::F64 mass;
    PS::F64vec pos_car;
    PS::F64vec pos;
    PS::F64ort boundary;
    template<class Tmom>
    void copyFromMoment(const Tmom & mom){
        this->mass     = mom.mass;
        this->pos      = mom.pos;
        this->pos_car        = mom.pos_car;
        this->boundary = mom.boundary;
    }
    void clear(){
        mass = 0.0;
        pos = 0.0;
        pos_car = 0.0;
        boundary.init();
    }
    PS::F64 getCharge() const {
        return mass;
    }
    PS::F64vec getPos() const {
        return pos;
    }
    PS::F64vec getPosCar() const {
        return pos_car;
    }
    MyMomentMonopole convertToMoment() const {
        return MyMomentMonopole(mass, pos, pos_car, boundary);
    }
};

class MyMomentQuadrupole{
public:
    PS::F64vec pos;
    PS::F64 mass;
    PS::F64mat quad;
    PS::F64vec pos_car;
    void init(){
        pos = 0.0;
        mass = 0.0;
        quad = 0.0;
        pos_car = 0.0;
    }
    MyMomentQuadrupole(){
        mass = 0.0;
        pos = 0.0;
        quad = 0.0;
        pos_car = 0.0;
    }
    MyMomentQuadrupole(const PS::F64 m, const PS::F64vec & p, const PS::F64mat & q, const PS::F64vec & p_car){
        mass = m;
        pos = p;
        quad = q;
        pos_car = p_car;
    }
    PS::F64vec getPos() const {
        return pos;
    }
    template<class Tepj>
    void accumulateAtLeaf(const Tepj & epj){
        mass += epj.getCharge();
        pos  += epj.getCharge() * epj.getPos();
        pos_car += epj.getCharge() * epj.getPosCar();
    }
    template<class Tepj>
    void accumulateAtLeaf2(const Tepj & epj){
        PS::F64 ctmp = epj.getCharge();
        PS::F64vec ptmp = epj.getPosCar() - this->pos_car;
        PS::F64 cx = ctmp * ptmp.x;
        PS::F64 cy = ctmp * ptmp.y;
        PS::F64 cz = ctmp * ptmp.z;
        this->quad.xx += cx * ptmp.x;
        this->quad.yy += cy * ptmp.y;
        this->quad.zz += cz * ptmp.z;
        this->quad.xy += cx * ptmp.y;
        this->quad.xz += cx * ptmp.z;
        this->quad.yz += cy * ptmp.z;
    }
    void set(){
        pos = pos / mass;
        pos_car = pos_car / mass;
    }
    void accumulate(const MyMomentQuadrupole & mom){
        mass += mom.mass;
        pos += mom.mass * mom.pos;
        pos_car += mom.mass * mom.pos_car;
    }
    void accumulate2(const MyMomentQuadrupole & mom){
        PS::F64 mtmp = mom.mass;
        PS::F64vec ptmp = mom.pos_car - this->pos_car;
        PS::F64 cx = mtmp * ptmp.x;
        PS::F64 cy = mtmp * ptmp.y;
        PS::F64 cz = mtmp * ptmp.z;
        this->quad.xx += cx * ptmp.x + mom.quad.xx;
        this->quad.yy += cy * ptmp.y + mom.quad.yy;
        this->quad.zz += cz * ptmp.z + mom.quad.zz;
        this->quad.xy += cx * ptmp.y + mom.quad.xy;
        this->quad.xz += cx * ptmp.z + mom.quad.xz;
        this->quad.yz += cy * ptmp.z + mom.quad.yz;
    }
    void dump(std::ostream & fout = std::cout) const {
        fout<<"mass= "<<mass<<std::endl;
        fout<<"pos= "<<pos<<std::endl;
        fout<<"quad= "<<quad<<std::endl;
    }
};

class MySPJQuadrupole{
public:
    PS::F64 mass;
    PS::F64vec pos;
    PS::F64vec pos_car;
    PS::F64mat quad;
    PS::F64 getCharge() const {
        return mass;
    }
    PS::F64vec getPos() const {
        return pos;
    }
    PS::F64vec getPosCar() const {
        return pos_car;
    }
    void copyFromMoment(const MyMomentQuadrupole & mom){
        mass = mom.mass;
        pos = mom.pos;
        quad = mom.quad;
        pos_car = mom.pos_car;
    }
    MyMomentQuadrupole convertToMoment() const {
        return MyMomentQuadrupole(mass, pos, quad, pos_car);
    }
    void clear(){
        mass = 0.0;
        pos = 0.0;
        quad = 0.0;
        pos_car = 0.0;
    }
};


#ifdef ENABLE_PHANTOM_GRAPE_X86


template <class TParticleJ>
void CalcGravity(const FPGrav * iptcl,
                 const PS::S32 ni,
                 const TParticleJ * jptcl,
                 const PS::S32 nj,
                 FPGrav * force) {
    const PS::S32 nipipe = ni;
    const PS::S32 njpipe = nj;
    PS::F64 (*xi)[3] = (PS::F64 (*)[3])malloc(sizeof(PS::F64) * nipipe * PS::DIMENSION);
    PS::F64 (*ai)[3] = (PS::F64 (*)[3])malloc(sizeof(PS::F64) * nipipe * PS::DIMENSION);
    PS::F64  *pi     = (PS::F64  *    )malloc(sizeof(PS::F64) * nipipe);
    PS::F64 (*xj)[3] = (PS::F64 (*)[3])malloc(sizeof(PS::F64) * njpipe * PS::DIMENSION);
    PS::F64  *mj     = (PS::F64  *    )malloc(sizeof(PS::F64) * njpipe);
    for(PS::S32 i = 0; i < ni; i++) {
        xi[i][0] = iptcl[i].getPos()[0];
        xi[i][1] = iptcl[i].getPos()[1];
        xi[i][2] = iptcl[i].getPos()[2];
        ai[i][0] = 0.0;
        ai[i][1] = 0.0;
        ai[i][2] = 0.0;
        pi[i]    = 0.0;
    }
    for(PS::S32 j = 0; j < nj; j++) {
        xj[j][0] = jptcl[j].getPos()[0];
        xj[j][1] = jptcl[j].getPos()[1];
        xj[j][2] = jptcl[j].getPos()[2];
        mj[j]    = jptcl[j].getCharge();
        xj[j][0] = jptcl[j].pos[0];
        xj[j][1] = jptcl[j].pos[1];
        xj[j][2] = jptcl[j].pos[2];
        mj[j]    = jptcl[j].mass;
    }
    PS::S32 devid = PS::Comm::getThreadNum();
    g5_set_xmjMC(devid, 0, nj, xj, mj);
    g5_set_nMC(devid, nj);
    g5_calculate_force_on_xMC(devid, xi, ai, pi, ni);
    for(PS::S32 i = 0; i < ni; i++) {
        force[i].acc[0] += ai[i][0];
        force[i].acc[1] += ai[i][1];
        force[i].acc[2] += ai[i][2];
        force[i].pot    -= pi[i];
    }
    free(xi);
    free(ai);
    free(pi);
    free(xj);
    free(mj);
}

#else



template <class TParticleJ>
void CalcForceEpold(const FPGrav * pi,
                 const PS::S32 ni,
                 const TParticleJ * pj,
                 const PS::S32 nj,
                 ForceGrav * force) {

#if 0
    std::cerr << "calcforceep  called with "<<ni << " "<< nj << " ips:\n";
    for(int i=0; i<ni; i++){
    	std::cerr << i << " "<< pi[i].pos_car<<
	    " "<< force[i].acc<< " " <<
	    force[i].pot<<"\n";
    }
    std::cerr << "jps:\n";
    
    for(int j=0; j<nj; j++){
    	std::cerr << j << " " << pj[j].mass << " " <<  pj[j].pos_car<< " "
		  <<  pj[j].pos<<
    	    "\n";
    }
#endif    

    PS::F64 eps2 = FPGrav::eps * FPGrav::eps;
    PS::F64 kappa = FPGrav::kappa;
    PS::F64 eta   = FPGrav::eta;
    PS::F64 xj[nj];
    PS::F64 yj[nj];
    PS::F64 zj[nj];
    PS::F64 vxj[nj];
    PS::F64 vyj[nj];
    PS::F64 vzj[nj];
    PS::F64 mj[nj];
    PS::S64 jid[nj];
    for(int j=0; j<nj; j++){
	xj[j] = pj[j].pos_car[0];
	yj[j] = pj[j].pos_car[1];
	zj[j] = pj[j].pos_car[2];
	vxj[j] = pj[j].vel[0];
	vyj[j] = pj[j].vel[1];
	vzj[j] = pj[j].vel[2];
	mj[j] = pj[j].getCharge();
	jid[j] = pj[j].id;
    }
	
    for(int i=0; i<ni; i++){
	const PS::F64vec xi = pi[i].pos_car;
	PS::F64vec ai = 0.0;
	PS::F64 poti = 0.0;
	PS::F64 dedtdumperi = 0.0;
	PS::F64 r_coll = FPGrav::rcoll*2;
	PS::F64  r_coll_inv = 1.0/r_coll;
	PS::F64 r_coll_sq = r_coll*r_coll;
	PS::F64 r_coll_third_inv = 1.0/(r_coll_sq*r_coll);
	PS::F64 r_coll_sq_inv = 1.0/(r_coll*r_coll);
	PS::F64 xix= pi[i].pos_car[0];
	PS::F64 xiy = pi[i].pos_car[1];
	PS::F64 xiz = pi[i].pos_car[2];
	PS::F64 aix = 0.0;
	PS::F64 aiy = 0.0;
	PS::F64 aiz = 0.0;
	PS::F64 vix = pi[i].vel[0];
	PS::F64 viy = pi[i].vel[1];
	PS::F64 viz = pi[i].vel[2];
	PS::F64 mi = pi[i].getCharge();
	PS::S64 iid = pi[i].id;
	
#pragma omp simd reduction(+:aix,aiy,aiz,poti,dedtdumperi)	    
	for(int j=0; j<nj; j++){
	    //	    PS::F64vec rij    = xi - pj[j].pos_car;
	    PS::F64 rijx    = xix - xj[j];
	    PS::F64 rijy    = xiy - yj[j];
	    PS::F64 rijz    = xiz - zj[j];
	    if(iid != jid[j]) {
		PS::F64 r2 = rijx*rijx +rijy*rijy +rijz*rijz +eps2;
		PS::F64 r2_inv  = 1.0/r2;
		PS::F64 r_inv  =  sqrt(r2_inv);
		PS::F64 pot = r_inv * mj[j];
#if 0	    
		fprintf(stderr, "i=%d j=%d r2 =%g r_coll_sq=%g\n",
			i, j, r2, r_coll_sq);
#endif	    
		if(r_coll_sq < r2){
		    aix     +=  -pot * r2_inv * rijx;
		    aiy     +=  -pot * r2_inv * rijy;
		    aiz     +=  -pot * r2_inv * rijz;
		    poti   +=  -pot;
		}else{
		    aix     +=  -rijx* mj[j]*r_coll_third_inv;;
		    aiy     +=  -rijy* mj[j]*r_coll_third_inv;;
		    aiz     +=  -rijz* mj[j]*r_coll_third_inv;;
		    poti   +=  -mj[j]*r_coll_inv;;
		    
		    PS::F64 m_r = mj[j]/ (mi+mj[j]);
		    PS::F64 r = sqrt(r2);
		    PS::F64 dr = r_coll-r ;
		    PS::F64 kappacoef = kappa * m_r * dr/r;
		    aix += kappacoef* rijx;
		    aiy += kappacoef* rijy;
		    aiz += kappacoef* rijz;
		    //poti += kappa*m_r*dr*dr * 0.5 -dr* pj[j].getCharge()*r_coll_sq_inv;
		    poti += kappacoef*dr*r * 0.5;
		    PS::F64vec vij = pi[i].vel - pj[j].vel;
		    PS::F64 vijx = vix - vxj[j];
		    PS::F64 vijy = viy - vyj[j];
		    PS::F64 vijz = viz - vzj[j];
		    PS::F64 rv = rijx*vijx + rijy*vijy + rijz*vijz;
		    PS::F64 etacoef = -eta * m_r * rv * r2_inv;
#if 0		
		    fprintf(stderr, " ri= %g %g %g vi= %g %g %g\n",
			    xix, xiy, xiz, vix, viy, viz);
		    fprintf(stderr, " dr = %g k force, pot =%g %g %g  %g\n",
			    dr, kappacoef* rijx, kappacoef* rijy,kappacoef* rijz,
			    kappa*dr*r * 0.5);
		    
		    fprintf(stderr, " eta=%g rv=%g  etacoef=%g f= %g %g %g dedt=%g\n",
			    eta, rv, etacoef,
			    etacoef * rijx, etacoef * rijy, etacoef * rijz,
			    -etacoef * rv);
#endif			
		    aix += etacoef * rijx;
		    aiy += etacoef * rijy;
		    aiz += etacoef * rijz;
		    dedtdumperi -= etacoef * rv*m_r;
		}
	    }
	}
	force[i].acc += PS::F64vec(aix, aiy, aiz);
	force[i].pot += poti;
	force[i].dedtdumper += dedtdumperi;
    }
    
}


#if 0
typedef PS::F64 SpFtype;    
    typedef  PS::F64 force_t;
typedef  PS::F64 sumforce_t;
    typedef  PS::S64 intforce_t;
#else
typedef  PS::F32 SpFtype;    
typedef  PS::F32 force_t;
typedef  PS::F32 sumforce_t;
typedef  PS::S32 intforce_t;
#endif

#ifdef USE_PIKG_KERNEL
//#if 0

struct PPEpi{
    PS::F32vec pos;
    PS::F32vec vel;
    PS::F32    mass;
    PS::S32    id;
};
struct PPEpj{
    PS::F32vec pos;
    PS::F32vec vel;
    PS::F32    mass;
    PS::S32    id;
};
struct PPForce{
    PS::F32vec acc;
    PS::F32    pot;
    PS::F32    dedtdumperi;
};

#include "epep.hpp"


template <class TParticleJ>
void CalcForceEp(const FPGrav * pi,
                 const PS::S32 ni,
                 const TParticleJ * pj,
                 const PS::S32 nj,
                 ForceGrav * force) {
    
#if 1
    std::cerr << "calcforceep  called with "<<ni << " "<< nj << " ips:\n";
    for(int i=0; i<ni; i++){
    	std::cerr << i << " "<< pi[i].pos_car<<
	    " "<< force[i].acc<< " " <<
	    force[i].pot<<"\n";
    }
    std::cerr << "jps:\n";
    
    for(int j=0; j<nj; j++){
    	std::cerr << j << " " << pj[j].mass << " " <<  pj[j].pos_car<< " "
		  <<  pj[j].pos<<
    	    "\n";
    }
#endif    

    
    PS::F32 eps2 = FPGrav::eps * FPGrav::eps+1e-30;
    
    PS::F32 kappa = FPGrav::kappa;
    PS::F32 eta   = FPGrav::eta;
    PS::F32 r_coll   = FPGrav::rcoll*2;
    auto r_coll_inv  = 1.0/r_coll;
    auto r_coll_sq  = r_coll*r_coll;
    
    
    PPEpi epi[ni];
    PPForce f[ni];
    for(int i=0;i<ni;i++){
	epi[i].pos = (PS::F32vec)(pi[i].pos_car - pi[0].pos_car);
	epi[i].vel = (PS::F32vec)(pi[i].vel - pi[0].vel);
	epi[i].mass = pi[i].mass;
	epi[i].id = pi[i].id;
	
	f[i].acc = 0.0;
	f[i].pot = 0.0;
	f[i].dedtdumperi = 0.0;
    }
    PPEpj epj[nj];
    for(int i=0;i<nj;i++){
	epj[i].pos = (PS::F32vec)(pj[i].pos_car - pi[0].pos_car);
	epj[i].vel = (PS::F32vec)(pj[i].vel - pi[0].vel);
	epj[i].mass = pj[i].mass;
	epj[i].id = pj[i].id;
	
    }
    epep_pikg(eps2, r_coll, r_coll_inv, r_coll_sq, kappa, eta)
	(epi,ni, epj,nj,f);
    for(int i=0;i<ni;i++){
#if 1
	fprintf(stderr, "Force[%d] = %e  %e  %e  %e %e\n",
		i, f[i].acc.x, f[i].acc.y, f[i].acc.z,
		f[i].pot, f[i].dedtdumperi);
#endif	    
	force[i].acc += f[i].acc;
	force[i].pot += f[i].pot;
	force[i].dedtdumper += f[i].dedtdumperi;
    }
}


#else

template <class TParticleJ>
void CalcForceEp(const FPGrav * pi,
                 const PS::S32 ni,
                 const TParticleJ * pj,
                 const PS::S32 nj,
                 ForceGrav * force) {

#if 0
    std::cerr << "calcforceep  called with "<<ni << " "<< nj << " ips:\n";
    for(int i=0; i<ni; i++){
    	std::cerr << i << " "<< pi[i].pos_car<<
	    " "<< force[i].acc<< " " <<
	    force[i].pot<<"\n";
    }
    std::cerr << "jps:\n";
    
    for(int j=0; j<nj; j++){
    	std::cerr << j << " " << pj[j].mass << " " <<  pj[j].pos_car<< " "
		  <<  pj[j].pos<<
    	    "\n";
    }
#endif    

    
    force_t eps2 = FPGrav::eps * FPGrav::eps+1e-30;
    force_t kappa = FPGrav::kappa;
    force_t eta   = FPGrav::eta;
    force_t xj[nj];
    force_t yj[nj];
    force_t zj[nj];
    force_t vxj[nj];
    force_t vyj[nj];
    force_t vzj[nj];
    force_t mj[nj];
    intforce_t jid[nj];
    force_t x0= pi[0].pos_car[0];
    force_t y0 = pi[0].pos_car[1];
    force_t z0 = pi[0].pos_car[2];
    force_t vx0 = pi[0].vel[0];
    force_t vy0 = pi[0].vel[1];
    force_t vz0 = pi[0].vel[2];
    
    for(int j=0; j<nj; j++){
	xj[j] = pj[j].pos_car[0] -x0;
	yj[j] = pj[j].pos_car[1]-y0;
	zj[j] = pj[j].pos_car[2]-z0;
	vxj[j] = pj[j].vel[0]-vx0;
	vyj[j] = pj[j].vel[1]-vy0;
	vzj[j] = pj[j].vel[2]-vz0;
	mj[j] = pj[j].getCharge();
	jid[j] = pj[j].id;
    }

    
    for(int i=0; i<ni; i++){
	//	const PS::F64vec xi = pi[i].pos_car;
	PS::F64vec ai = 0.0;
	sumforce_t poti = 0.0;
	force_t dedtdumperi = 0.0;
 	force_t r_coll = FPGrav::rcoll*2;
	force_t  r_coll_inv = 1.0/r_coll;
	force_t r_coll_sq = r_coll*r_coll;
	force_t r_coll_third_inv = 1.0/(r_coll_sq*r_coll);
	force_t r_coll_sq_inv = 1.0/(r_coll*r_coll);
	force_t xix= pi[i].pos_car[0] -x0;
	force_t xiy = pi[i].pos_car[1] -y0;
	force_t xiz = pi[i].pos_car[2] -z0;
	sumforce_t aix = 0.0;
	sumforce_t aiy = 0.0;
	sumforce_t aiz = 0.0;
	force_t vix = pi[i].vel[0] -vx0;
	force_t viy = pi[i].vel[1] -vy0;
	force_t viz = pi[i].vel[2] -vz0;
	force_t mi = pi[i].getCharge();
	intforce_t iid = pi[i].id;
	
#pragma omp simd reduction(+:aix,aiy,aiz,poti,dedtdumperi)	    
	for(int j=0; j<nj; j++){
	    //	    PS::F64vec rij    = xi - pj[j].pos_car;
	    force_t rijx    = xix - xj[j];
	    force_t rijy    = xiy - yj[j];
	    force_t rijz    = xiz - zj[j];
	    force_t mjlocal = mj[j];
	    force_t m_r = mjlocal/ (mi+mjlocal);
	    force_t r2 = rijx*rijx +rijy*rijy +rijz*rijz +eps2;
	    force_t r2_inv;
	    force_t r_inv;
	    force_t r = sqrt(r2);
	    force_t r_inv_true = 1.0/r;
	    auto r2_inv_true =  r_inv_true*r_inv_true;
	    force_t dr = r_coll-r ;
	    if(r_coll_sq < r2){
		r_inv = r_inv_true;
		r2_inv = r2_inv_true;
		dr = 0.0;
		m_r = 0.0;
	    }else{
		//	fprintf(stderr, "in collision: %d %d\n", iid, jid[j]);
		r_inv = r_coll_inv;
		r2_inv = r_inv*r_inv;
		if(iid == jid[j]) {
		    m_r=0;
		    mjlocal=0;
		}
		
	    }
	    force_t pot = r_inv * mjlocal;
#if 0	    
	    fprintf(stderr, "i=%d j=%d r2 =%g r_coll_sq=%g\n",
		    i, j, r2, r_coll_sq);
#endif	    
	    aix     +=  -pot * r2_inv * rijx;
	    aiy     +=  -pot * r2_inv * rijy;
	    aiz     +=  -pot * r2_inv * rijz;
	    poti   +=  -pot;
	    force_t kappacoef = kappa * m_r * dr * r_inv_true;
	    //	    fprintf(stderr, "r2=%e r_coll=%e  kappac=%e pot=%e %e %e\n",
	    //		    r2, r_coll, kappacoef, pot, r2_inv, pot * r2_inv * rijx);
	    aix += kappacoef* rijx;
	    aiy += kappacoef* rijy;
	    aiz += kappacoef* rijz;
	    //poti += kappa*m_r*dr*dr * 0.5 -dr* pj[j].getCharge()*r_coll_sq_inv;
	    poti += kappacoef*dr*r * 0.5;
	    //	    PS::F64vec vij = pi[i].vel - pj[j].vel;
	    force_t vijx = vix - vxj[j];
	    force_t vijy = viy - vyj[j];
	    force_t vijz = viz - vzj[j];
	    force_t rv = rijx*vijx + rijy*vijy + rijz*vijz;
	    if(r_coll_sq < r2) rv=0.0;
	    
	    force_t etacoef = -eta * m_r * rv * r2_inv_true;
#if 0		
	    fprintf(stderr, " ri= %g %g %g vi= %g %g %g\n",
		    xix, xiy, xiz, vix, viy, viz);
	    fprintf(stderr, " dr = %g k force, pot =%g %g %g  %g\n",
		    dr, kappacoef* rijx, kappacoef* rijy,kappacoef* rijz,
		    kappa*dr*r * 0.5);
	    
	    fprintf(stderr, " eta=%g rv=%g  etacoef=%g f= %g %g %g dedt=%g\n",
		    eta, rv, etacoef,
		    etacoef * rijx, etacoef * rijy, etacoef * rijz,
		    -etacoef * rv);
#endif			
	    aix += etacoef * rijx;
	    aiy += etacoef * rijy;
	    aiz += etacoef * rijz;

	    dedtdumperi -= etacoef * rv*m_r;
	    //	    fprintf(stderr, "%d %d %e %e\n",
	    //		    iid, jid[j], aix, aiy);
	}
#if 0
	fprintf(stderr, "Force[%d] = %e  %e  %e  %e %e\n",
		i, aix, aiy, aiz, poti, dedtdumperi);
#endif	    
	force[i].acc += PS::F64vec(aix, aiy, aiz);
	force[i].pot += poti;
	force[i].dedtdumper += dedtdumperi;
    }
    
}

#endif


template <class TParticleJ>
void CalcGravitySp(const FPGrav * ep_i,
		   const PS::S32 n_ip,
		   const TParticleJ * ep_j,
		   const PS::S32 n_jp,
		   FPGrav * force) {
    //    std::cerr << "calcgravitysp  called with "<<n_ip << " "<< n_jp << "\n";
    PS::F64 eps2 = FPGrav::eps * FPGrav::eps;
    for(PS::S32 i = 0; i < n_ip; i++){
        PS::F64vec xi = ep_i[i].pos;
        PS::F64vec ai = 0.0;
        PS::F64 poti = 0.0;
        for(PS::S32 j = 0; j < n_jp; j++){
	    PS::F64vec rij    = xi - ep_j[j].pos;
	    PS::F64    r3_inv = rij * rij + eps2;
	    PS::F64    r_inv  = 1.0/sqrt(r3_inv);
	    r3_inv  = r_inv * r_inv;
	    r_inv  *= ep_j[j].getCharge();
	    r3_inv *= r_inv;
	    ai     -= r3_inv * rij;
	    poti   -= r_inv;
	}
	force[i].acc += ai;
	force[i].pot += poti;
    }
}

#endif

template<typename Tpi, typename Tpj, typename Tforce>
struct CalcForceSpMono{
    void operator ()(const Tpi * pi,
                     const PS::S32 ni,
                     const Tpj * pj,
                     const PS::S32 nj,
                     Tforce * force){
	//std::cerr << "calcuforcespmono called with "<<ni << " "<< nj << "\n";
        const auto eps2 = FPGrav::eps*FPGrav::eps;
	PS::F64vec xj[nj];

#if 0	
	std::cerr << "calcforcesomono  called with "<<ni << " "<< nj << " ips:\n";
	for(int i=0; i<ni; i++){
	    std::cerr << i << " "<< pi[i].getPosCar()<< " " <<
		pi[i].getPos() << " "
		      << force[i].acc<< " " <<
		force[i].pot<<"\n";"\n";
	}
	std::cerr << "jps:\n";
	
	for(int j=0; j<nj; j++){
	    std::cerr << j << " " << pj[j].getCharge() << " " <<  pj[j].getPosCar()<< " "
		      <<  pj[j].getPos()<<    	    "\n";
	}
#endif	
	for(auto j=0; j<nj; j++){
	    xj[j] = pj[j].getPosCar();
	}
	for(auto i=0; i<ni; i++){
	    const auto xi = pi[i].getPosCar();
	    PS::F64vec ai = 0.0;
	    PS::F64 poti  = 0.0;
	    for(auto j=0; j<nj; j++){
		//PS::F64vec rij    = xi - pj[j].getPosCar();
		const auto rij    = xi - xj[j];
		auto r3_inv = rij * rij + eps2;
		auto r_inv  = 1.0/sqrt(r3_inv);
		r3_inv  = r_inv * r_inv;
		r_inv  *= pj[j].getCharge();
		r3_inv *= r_inv;
		ai     -= r3_inv * rij;
		poti   -= r_inv;
		//std::cerr << i<< " " << j << " " << r_inv << " " << poti <<"\n";
	    }
	    force[i].acc += ai;
	    force[i].pot += poti;
	}
#if 0	
	std::cerr << "calcforcesomono  results "<<ni << " "<< nj << " ips:\n";
	for(int i=0; i<ni; i++){
	    std::cerr << i << " "<< force[i].acc<< " " <<
		force[i].pot<<"\n";
	}
#endif	
	
    }
};

template<typename Tpi, typename Tpj, typename Tforce>
struct CalcForceSpQuadold{
    void operator ()(const Tpi * pi,
                     const PS::S32 ni,
                     const Tpj * pj,
                     const PS::S32 nj,
                     Tforce * force){
        const auto eps2 = FPGrav::eps*FPGrav::eps;
	PS::F64vec xj[nj];
	for(auto j=0; j<nj; j++){
	    xj[j] = pj[j].getPosCar();
	}
        for(auto i=0; i<ni; i++){
            PS::F64vec xi = pi[i].getPosCar();
            PS::F64vec ai = 0.0;
            PS::F64 poti = 0.0;
            for(auto j=0; j<nj; j++){
                PS::F64 mj = pj[j].getCharge();
                //PS::F64vec xj= pj[j].getPosCar();
                //PS::F64vec rij= xi - xj;
		PS::F64vec rij = xi - xj[j];
                PS::F64 r2 = rij * rij + eps2;
                PS::F64mat qj = pj[j].quad;
                PS::F64 tr = qj.getTrace();
                PS::F64vec qr( (qj.xx*rij.x + qj.xy*rij.y + qj.xz*rij.z),
                               (qj.yy*rij.y + qj.yz*rij.z + qj.xy*rij.x),
                               (qj.zz*rij.z + qj.xz*rij.x + qj.yz*rij.y) );
                PS::F64 qrr = qr * rij;
                PS::F64 r_inv = 1.0f/sqrt(r2);
                PS::F64 r2_inv = r_inv * r_inv;
                PS::F64 r3_inv = r2_inv * r_inv;
                PS::F64 r5_inv = r2_inv * r3_inv * 1.5;
                PS::F64 qrr_r5 = r5_inv * qrr;
                PS::F64 qrr_r7 = r2_inv * qrr_r5;
                PS::F64 A = mj*r3_inv - tr*r5_inv + 5*qrr_r7;
                PS::F64 B = -2.0*r5_inv;
                ai -= A*rij + B*qr;
                poti -= mj*r_inv - 0.5*tr*r3_inv + qrr_r5;
            }
            force[i].acc += ai;
            force[i].pot += poti;
        }
    }
};

template<typename Tpi, typename Tpj, typename Tforce>
struct CalcForceSpQuadold2{
    void operator ()(const Tpi * pi,
                     const PS::S32 ni,
                     const Tpj * pj,
                     const PS::S32 nj,
                     Tforce * force){
#if 0	
	std::cerr << "calcforcesspquad  called with "<<ni << " "<< nj << " ips:\n";
	for(int i=0; i<ni; i++){
	    std::cerr << i << " "<< pi[i].getPosCar()<< " " <<
		pi[i].getPos() << " "
		      << force[i].acc<< " " <<
		force[i].pot<<"\n";"\n";
	}
	std::cerr << "jps:\n";
	
	for(int j=0; j<nj; j++){
	    std::cerr << j << " " << pj[j].getCharge() << " " <<  pj[j].getPosCar()<< " "
		      <<  pj[j].getPos()<<    	    "\n";
	}
#endif	
	const auto eps2 = FPGrav::eps*FPGrav::eps;
	PS::F64 xjx[nj];
	PS::F64 xjy[nj];
	PS::F64 xjz[nj];
	PS::F64 mj[nj];
	PS::F64 qjxx[nj];
	PS::F64 qjxy[nj];
	PS::F64 qjxz[nj];
	PS::F64 qjyy[nj];
	PS::F64 qjyz[nj];
	PS::F64 qjzz[nj];
	PS::F64 tr[nj];
	for(auto j=0; j<nj; j++){
	    xjx[j] = pj[j].getPosCar()[0];
	    xjy[j] = pj[j].getPosCar()[1];
	    xjz[j] = pj[j].getPosCar()[2];
	    mj[j] = pj[j].getCharge();
	    qjxx[j] = pj[j].quad.xx;
	    qjxy[j] = pj[j].quad.xy;
	    qjxz[j] = pj[j].quad.xz;
	    qjyy[j] = pj[j].quad.yy;
	    qjyz[j] = pj[j].quad.yz;
	    qjzz[j] = pj[j].quad.zz;
	    tr[j] = pj[j].quad.getTrace();
	}
	for(auto i=0; i<ni; i++){
	    PS::F64 xix = pi[i].getPosCar()[0];
	    PS::F64 xiy = pi[i].getPosCar()[1];
	    PS::F64 xiz = pi[i].getPosCar()[2];

	    PS::F64 aix = 0.0;
	    PS::F64 aiy = 0.0;
	    PS::F64 aiz = 0.0;
	    PS::F64 poti = 0.0;
#pragma clan loop unroll_count(6)
#pragma omp simd reduction(+:aix,aiy,aiz,poti)	    
	    for(auto j=0; j<nj; j++){
		//		PS::F64vec rij = xi - xj[j];
		PS::F64 rijx = xix - xjx[j];
		PS::F64 rijy = xiy - xjy[j];
		PS::F64 rijz = xiz - xjz[j];
		PS::F64 r2 = rijx * rijx +rijy * rijy +rijz * rijz +eps2;
		PS::F64 qrx=qjxx[j]*rijx + qjxy[j]*rijy + qjxz[j]*rijz;
		PS::F64 qry=qjyy[j]*rijy + qjyz[j]*rijz + qjxy[j]*rijx;
		PS::F64 qrz=qjzz[j]*rijz + qjxz[j]*rijx + qjyz[j]*rijy;
		PS::F64 qrr = qrx*rijx +qry*rijy +qrz*rijz;
		PS::F64 r_inv = 1.0f/sqrt(r2);
		PS::F64 r2_inv = r_inv * r_inv;
		PS::F64 r3_inv = r2_inv * r_inv;
		PS::F64 r5_inv = r2_inv * r3_inv * 1.5;
		PS::F64 qrr_r5 = r5_inv * qrr;
		PS::F64 qrr_r7 = r2_inv * qrr_r5;
		PS::F64 A = mj[j]*r3_inv - tr[j]*r5_inv + 5*qrr_r7;
		PS::F64 B = -2.0*r5_inv;
		aix += A*rijx + B*qrx;
		aiy += A*rijy + B*qry;
		aiz += A*rijz + B*qrz;
		poti += mj[j]*r_inv - 0.5*tr[j]*r3_inv + qrr_r5;
	    }
	    force[i].acc -= PS::F64vec(aix, aiy, aiz);
	    force[i].pot -= poti;
	}
    }
};


#ifdef USE_PIKG_KERNEL
struct Epi{
    PS::F32vec pos;
};
struct Epj{
    PS::F32vec pos;
    PS::F32    mass;
    PS::F32    qxx;
    PS::F32    qyy;
    PS::F32    qzz;
    PS::F32    qxy;
    PS::F32    qyz;
    PS::F32    qxz;
    PS::F32    qtr;
};
struct Force{
  PS::F32vec acc;
  PS::F32    pot;
};

#include "epsp.hpp"



template<typename Tpi, typename Tpj, typename Tforce>
struct CalcForceSpQuad{
    void operator ()(const Tpi * pi,
                     const PS::S32 ni,
                     const Tpj * pj,
                     const PS::S32 nj,
                     Tforce * force){
#if 0
	std::cerr << "calcforcesspquad  called with "<<ni << " "<< nj << " ips:\n";
	for(int i=0; i<ni; i++){
	    std::cerr << i << " "<< pi[i].getPosCar()<< " " <<
		pi[i].getPos() << " "
		      << force[i].acc<< " " <<
		force[i].pot<<"\n";"\n";
	}
	std::cerr << "jps:\n";
	
	for(int j=0; j<nj; j++){
	    std::cerr << j << " " << pj[j].getCharge() << " " <<  pj[j].getPosCar()<< " "
		      <<  pj[j].getPos()<<    	    "\n";
	}
#endif	
	const auto eps2 = FPGrav::eps*FPGrav::eps;
	Epi epi[ni];
	Force f[ni];
	for(int i=0;i<ni;i++){
	    epi[i].pos = (PS::F32vec)(pi[i].pos_car - pi[0].pos_car);
	    
	    f[i].acc = 0.0;
	    f[i].pot = 0.0;
	}
	Epj epj[nj];
	for(int i=0;i<nj;i++){
	    epj[i].pos = (PS::F32vec)(pj[i].pos_car - pi[0].pos_car);
	    epj[i].mass = pj[i].mass;
	    epj[i].qxx = pj[i].quad.xx;
	    epj[i].qxy = pj[i].quad.xy;
	    epj[i].qxz = pj[i].quad.xz;
	    epj[i].qyy = pj[i].quad.yy;
	    epj[i].qyz = pj[i].quad.yz;
	    epj[i].qzz = pj[i].quad.zz;
	    epj[i].qtr = pj[i].quad.getTrace();
	    
	}
	epsp_pikg((PIKG::F32)eps2)(epi,ni, epj,nj,f);
	for(int i=0;i<ni;i++){
#if 0
	    fprintf(stderr, "Force[%d] = %e  %e  %e  %e\n",
		    i, f[i].acc.x, f[i].acc.y, f[i].acc.z,
		    f[i].pot);
#endif	    
	    force[i].acc -= f[i].acc;
	    force[i].pot -= f[i].pot;
	}
    }
};

#else // non_pikg

template<typename Tpi, typename Tpj, typename Tforce>
struct CalcForceSpQuad{
    void operator ()(const Tpi * pi,
                     const PS::S32 ni,
                     const Tpj * pj,
                     const PS::S32 nj,
                     Tforce * force){
#if 0	
	std::cerr << "calcforcesspquad  called with "<<ni << " "<< nj << " ips:\n";
	for(int i=0; i<ni; i++){
	    std::cerr << i << " "<< pi[i].getPosCar()<< " " <<
		pi[i].getPos() << " "
		      << force[i].acc<< " " <<
		force[i].pot<<"\n";"\n";
	}
	std::cerr << "jps:\n";
	
	for(int j=0; j<nj; j++){
	    std::cerr << j << " " << pj[j].getCharge() << " " <<  pj[j].getPosCar()<< " "
		      <<  pj[j].getPos()<<    	    "\n";
	}
#endif	
	const auto eps2 = FPGrav::eps*FPGrav::eps;
	SpFtype xjx[nj];
	SpFtype xjy[nj];
	SpFtype xjz[nj];
	SpFtype mj[nj];
	SpFtype qjxx[nj];
	SpFtype qjxy[nj];
	SpFtype qjxz[nj];
	SpFtype qjyy[nj];
	SpFtype qjyz[nj];
	SpFtype qjzz[nj];
	SpFtype tr[nj];
	PS::F64 xix0 = pi[0].getPosCar()[0];
	PS::F64 xiy0 = pi[0].getPosCar()[1];
	PS::F64 xiz0 = pi[0].getPosCar()[2];
	for(auto j=0; j<nj; j++){
	    xjx[j] = pj[j].getPosCar()[0]-xix0;
	    xjy[j] = pj[j].getPosCar()[1]-xiy0;
	    xjz[j] = pj[j].getPosCar()[2]-xiz0;
	    mj[j] = pj[j].getCharge();
	    qjxx[j] = pj[j].quad.xx;
	    qjxy[j] = pj[j].quad.xy;
	    qjxz[j] = pj[j].quad.xz;
	    qjyy[j] = pj[j].quad.yy;
	    qjyz[j] = pj[j].quad.yz;
	    qjzz[j] = pj[j].quad.zz;
	    tr[j] = pj[j].quad.getTrace();
	}
	for(auto i=0; i<ni; i++){
	    SpFtype xix = pi[i].getPosCar()[0]-xix0;
	    SpFtype xiy = pi[i].getPosCar()[1]-xiy0;
	    SpFtype xiz = pi[i].getPosCar()[2]-xiz0;

	    SpFtype aix = 0.0;
	    SpFtype aiy = 0.0;
	    SpFtype aiz = 0.0;
	    SpFtype poti = 0.0;
#pragma clan loop unroll_count(6)
#pragma omp simd reduction(+:aix,aiy,aiz,poti)	    
	    for(auto j=0; j<nj; j++){
		//		PS::F64vec rij = xi - xj[j];
		SpFtype rijx = xix - xjx[j];
		SpFtype rijy = xiy - xjy[j];
		SpFtype rijz = xiz - xjz[j];
		SpFtype r2 = rijx * rijx +rijy * rijy +rijz * rijz +eps2;
		SpFtype qrx=qjxx[j]*rijx + qjxy[j]*rijy + qjxz[j]*rijz;
		SpFtype qry=qjyy[j]*rijy + qjyz[j]*rijz + qjxy[j]*rijx;
		SpFtype qrz=qjzz[j]*rijz + qjxz[j]*rijx + qjyz[j]*rijy;
		SpFtype qrr = qrx*rijx +qry*rijy +qrz*rijz;
		SpFtype r_inv = 1.0f/sqrt(r2);
		SpFtype r2_inv = r_inv * r_inv;
		SpFtype r3_inv = r2_inv * r_inv;
		SpFtype r5_inv = r2_inv * r3_inv * 1.5;
		SpFtype qrr_r5 = r5_inv * qrr;
		SpFtype qrr_r7 = r2_inv * qrr_r5;
		SpFtype A = mj[j]*r3_inv - tr[j]*r5_inv + 5*qrr_r7;
		SpFtype B = -2.0*r5_inv;
		aix += A*rijx + B*qrx;
		aiy += A*rijy + B*qry;
		aiz += A*rijz + B*qrz;
		poti += mj[j]*r_inv - 0.5*tr[j]*r3_inv + qrr_r5;
	    }
#if 0
	    fprintf(stderr, "Force[%d] = %e  %e  %e  %e\n",
		    i, aix, aiy, aiz, poti);
#endif	    
	    force[i].acc -= PS::F64vec(aix, aiy, aiz);
	    force[i].pot -= poti;
	}
    }
};

#endif

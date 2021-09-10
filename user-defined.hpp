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

class FPGrav{
public:
    PS::S64    id;
    PS::F64    mass;
    PS::F64vec pos;
    PS::F64vec pos_car;
    PS::F64vec vel;
    PS::F64vec acc;
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
	const auto pos_r   = sqrt(pos_car.x*pos_car.x + pos_car.y*pos_car.y);
	auto pos_phi = atan2(pos_car.y, pos_car.x);
	if (pos_phi >= M_PI) pos_phi -= SMALL;
	if (pos_phi <= -M_PI) pos_phi += SMALL;
	pos =  PS::F64vec(pos_phi, pos_r, pos_car.z);
    }
    void ctod()
    {
	const auto cth = cos(pos.x);
	const auto sth = sin(pos.x);
	const auto r = pos.y;
	const auto pos_x = r*cth;
	const auto pos_y = r*sth;
	pos_car= PS::F64vec(pos_x, pos_y, pos.z);
    }
    
	
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

    void copyFromForce(const FPGrav & force) {
        acc = force.acc;
        pot = force.pot;
    }

    void clear() {
        acc = 0.0;
        pot = 0.0;
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
void CalcForceEpworking(const FPGrav * ep_i,
                 const PS::S32 n_ip,
                 const TParticleJ * ep_j,
                 const PS::S32 n_jp,
                 FPGrav * force) {
    PS::F64 eps2 = FPGrav::eps * FPGrav::eps;
    for(PS::S32 i = 0; i < n_ip; i++){
        PS::F64vec xi = ep_i[i].getPos();
        PS::F64vec ai = 0.0;
        PS::F64 poti = 0.0;
        for(PS::S32 j = 0; j < n_jp; j++){
	    PS::F64vec rij    = xi - ep_j[j].getPos();
	    PS::F64    r2 = rij * rij + eps2;
	    if (r2 > 0){
		PS::F64    r2_inv  = 1.0/r2;
		PS::F64    r_inv  = sqrt(r2_inv);
		PS::F64 mr_inv = r_inv* ep_j[j].getCharge();
		PS::F64  mr3_inv = mr_inv*r2_inv;
		ai     -= mr3_inv * rij;
		poti   -= mr_inv;
	    }
	}
	force[i].acc += ai;
	force[i].pot += poti;
    }
}
template <class TParticleJ>
void CalcForceEporiginal(const FPGrav * ep_i,
                 const PS::S32 n_ip,
                 const TParticleJ * ep_j,
                 const PS::S32 n_jp,
                 FPGrav * force) {
    PS::F64 eps2 = FPGrav::eps * FPGrav::eps;
    for(PS::S32 i = 0; i < n_ip; i++){
        PS::F64vec xi = ep_i[i].getPos();
        PS::F64vec ai = 0.0;
        PS::F64 poti = 0.0;
        for(PS::S32 j = 0; j < n_jp; j++){
            PS::F64vec rij    = xi - ep_j[j].getPos();
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


template <class TParticleJ>
void CalcForceEp(const FPGrav * pi,
                 const PS::S32 ni,
                 const TParticleJ * pj,
                 const PS::S32 nj,
                 FPGrav * force) {
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
    for(int i=0; i<ni; i++){
	const PS::F64vec xi = pi[i].pos_car;
	PS::F64vec ai = 0.0;
	PS::F64vec ai_dash = 0.0;
	PS::F64 poti = 0.0;
	PS::F64 r_coll = FPGrav::rcoll;
	PS::F64  r_coll_inv = 1.0/r_coll;
	PS::F64 r_coll_sq = r_coll*r_coll;
	PS::F64 r_coll_third_inv = 1.0/(r_coll_sq*r_coll);
	PS::F64 r_coll_sq_inv = 1.0/(r_coll*r_coll);
	
	for(int j=0; j<nj; j++){
	    PS::F64vec rij    = xi - pj[j].pos_car;
	     if(pi[i].id == pj[j].id) continue;
	    PS::F64 r2 = rij * rij + eps2;
	    PS::F64 r2_inv  = 1.0/r2;
	    PS::F64 r_inv  =  sqrt(r2_inv);
	    PS::F64 pot = r_inv * pj[j].getCharge();
	    if(r_coll_sq < r2){
		ai     -= pot * r2_inv * rij;
		poti   -=  pot;
		//std::cerr << i<< " " << j << " " << pot << " " << poti << "\n";
	    }else{
		ai     -=  rij* pj[j].getCharge()*r_coll_third_inv;;
		poti   -=  pj[j].getCharge()*r_coll_inv;;
		
		//		fprintf(stderr, "in collision %lld %lld %g\n", pi[i].id,pj[j].id, r2);
		PS::F64 m_r = pj[j].mass / (pi[i].mass+pj[j].mass);
		PS::F64 r = sqrt(r2);
		PS::F64 dr = r_coll-r ;
		PS::F64vec a_kappa = kappa * m_r * dr/r * rij;
		ai += a_kappa;
		poti += kappa*m_r*dr*dr * 0.5 -dr* pj[j].getCharge()*r_coll_sq_inv;
		PS::F64vec vij = pi[i].vel - pj[j].vel;
		PS::F64 rv = rij*vij;
		PS::F64vec a_eta = -eta * m_r * rv * r2_inv * rij;
		// fprintf(stderr, "dr x,y,z, fx, fy, fz = %g %g %g %g %g %g %g %g %g %g\n",
		// 	dr, pi[i].pos.x,pi[i].pos.y, pi[i].pos.z,
		// 	a_kappa.x, a_kappa.y, a_kappa.z,
		// 	a_eta.x, a_eta.y, a_eta.z);
		
		ai_dash += a_eta;
		ai += a_eta;
	    }
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
struct CalcForceSpQuad{
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


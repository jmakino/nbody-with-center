#pragma once
class FileHeader{
public:
    PS::S64 n_body;
    PS::F64 time;
    PS::S32 readAscii(FILE * fp) {
        fscanf(fp, "%lf\n", &time);
        fscanf(fp, "%lld\n", &n_body);
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

    PS::F64 getCharge() const {
        return mass;
    }

    void copyFromFP(const FPGrav & fp){ 
        id = fp.id;
        mass = fp.mass;
        pos  = fp.pos;
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
                this->pos.x, this->pos.y, this->pos.z,
                this->vel.x, this->vel.y, this->vel.z);
    }

    void readAscii(FILE* fp) {
        fscanf(fp, "%lld\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", 
               &this->id, &this->mass,
               &this->pos.x, &this->pos.y, &this->pos.z,
               &this->vel.x, &this->vel.y, &this->vel.z);
        }

    void readBinary(FILE* fp) {       }

    void add_central_gravity()
    {
	PS::F64 r2=pos*pos;
	acc -= FPGrav::mass_center* pos/(r2*sqrt(r2));
	pot_center = -FPGrav::mass_center/(sqrt(r2));
	//	fprintf(stderr,"pos = %g %g %g r2= %g pot=%g pot_center= %g\n",
	//	pos.x, pos.y, pos.z, r2, pot, pot_center);
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
    PS::F64 eps2 = FPGrav::eps * FPGrav::eps;
    const auto kappa = FPGrav::kappa;
    const auto eta   = FPGrav::eta;
    for(auto i=0; i<ni; i++){
	const PS::F64vec xi = pi[i].getPos();
	PS::F64vec ai = 0.0;
	PS::F64vec ai_dash = 0.0;
	PS::F64 poti = 0.0;
	const auto r_coll = FPGrav::rcoll;
	PS::F64  r_coll_inv = 1.0/r_coll;
	const auto r_coll_sq = r_coll*r_coll;
	PS::F64 r_coll_third_inv = 1.0/(r_coll_sq*r_coll);
	PS::F64 r_coll_sq_inv = 1.0/(r_coll*r_coll);
	
	for(auto j=0; j<nj; j++){
	    PS::F64vec rij    = xi - pj[j].getPos();
	     if(pi[i].id == pj[j].id) continue;
	    PS::F64 r2 = rij * rij + eps2;
	    PS::F64 r2_inv  = 1.0/r2;
	    PS::F64 r_inv  =  sqrt(r2_inv);
	    PS::F64 pot = r_inv * pj[j].getCharge();
	    if(r_coll_sq < r2){
		ai     -= pot * r2_inv * rij;
		poti   -=  pot;
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
		PS::F64vec a_eta = eta * m_r * rv * r2_inv * rij;
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
}

#if 0
template<class Tpj>
CalcForceEp(const FPGrav * pi,
	    const PS::S32 ni,
	    const Tpj * pj,
	    const PS::S32 nj,
	    FPGrav * force){
    const auto eps2 = FPGrav::eps*FPGrav::eps;
    const auto kappa = FPGrav::kappa;
    const auto eta   = FPGrav::eta;
    for(auto i=0; i<ni; i++){
	const PS::F64vec xi = pi[i].getPos();
	PS::F64vec ai = 0.0;
	PS::F64vec ai_dash = 0.0;
	PS::F64 poti = 0.0;
	
	for(auto j=0; j<nj; j++){
	    const auto r_coll = FPGrav::rcoll;
	    const auto r_coll_sq = r_coll*r_coll;
	    PS::F64vec rij    = xi - pj[j].getPos();
	    if(pi[i].id == pj[j].id) continue;
	    PS::F64 r2 = rij * rij + eps2;
	    PS::F64 r_inv  = 1.0/sqrt(r2);
	    PS::F64 r2_inv  = r_inv * r_inv;
	    PS::F64 pot = r_inv * pj[j].getCharge();
	    ai     -= pot * r2_inv * rij;
	    poti   -= 0.5 * pot;
	    if(r_coll_sq > r2){
		PS::F64 m_r = pj[j].mass / (pi[i].mass+pj[j].mass);
		PS::F64 r = sqrt(r2);
		PS::F64 dr = r_coll-r ;
		ai += kappa * m_r * dr/r * rij;
		poti += 0.5*kappa*m_r*dr*dr * 0.5;
		PS::F64vec vij = pi[i].vel - pj[j].vel;
		PS::F64 rv = rij*vij;
		PS::F64vec a_eta = eta * m_r * rv * r2_inv * rij;
		ai_dash += a_eta;
		ai += a_eta;
	    }
	}
	force[i].acc += ai;
	force[i].pot += poti;
    }
}

#endif

template <class TParticleJ>
void CalcGravitySp(const FPGrav * ep_i,
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

#endif

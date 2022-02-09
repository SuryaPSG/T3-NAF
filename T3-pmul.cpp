#include <iostream>
#include <string>
#include <fstream>
#include <gmp.h>
#include <gmpxx.h>
#include <complex>
#include <NTL/GF2X.h>
#include <NTL/ZZ.h>
#include <vector>
#include <bitset>
#include <algorithm>
#include <unordered_map>
#include <chrono>

using namespace std;
using namespace NTL;
using namespace std::chrono;

//Function to convert hex string to GF2X

void hextoGF2X(GF2X& a,string inp){
	int len=inp.length(),t;
	for(int i=len-1;i>=0;i--){
		if(inp[i]>='0' && inp[i]<='9') t=inp[i]-'0';
		else if(inp[i]>='a' && inp[i]<='f') t=inp[i]-'a'+10;
		else if(inp[i]>='A' && inp[i]<='F') t=inp[i]-'A'+10;
		bitset<4> b(t);
		for(int j=0;j<4;j++){
			if(b[j]==1)
			SetCoeff(a,4*(len-i-1)+j,1);
			else
			SetCoeff(a,4*(len-i-1)+j,0);
		}
	}
}

class LDPoint;

class Point{

//Class for holding the points on a curve

public:
	GF2X x;
	GF2X y;
	Point(){
		memset(this, 0, sizeof(Point));
	}
	void toLDPoint(LDPoint& P);
	void negate(Point& R,const GF2X mt){
		R.x=(this->x)%mt;
		R.y=(this->x+this->y)%mt;
		
	}
	void print(){
		GF2X::HexOutput=1;
		cout<<"x : "<<this->x<<endl;
		cout<<"y : "<<this->y<<endl;
		GF2X::HexOutput=0;
	}
};

class LDPoint{

//Point in Lopez-Dahab co-ordinate system (X,Y,Z)

public :
	GF2X X;
	GF2X Y;
	GF2X Z;
	LDPoint(){
		memset(this, 0, sizeof(LDPoint));
	}
	void toPoint(Point& P,const GF2X mt){
		if(!IsZero(this->Z)){
			GF2X temp;
			temp=InvMod(this->Z,mt);
			P.x=(this->X*temp)%mt;
			P.y=(this->Y*(temp*temp)%mt)%mt;
		}
	}
	void square(const GF2X& mt){
		this->X=SqrMod(X,mt);
		this->Y=SqrMod(Y,mt);
		this->Z=SqrMod(Z,mt);
	}
	void set(const LDPoint& P){
		this->X=P.X;
		this->Y=P.Y;
		this->Z=P.Z;
	}
	void print(){
		GF2X::HexOutput=1;
		cout<<"X : "<<this->X<<endl;
		cout<<"Y : "<<this->Y<<endl;
		cout<<"Z : "<<this->Z<<endl;
		GF2X::HexOutput=0;
	}
};

void Point::toLDPoint(LDPoint& P){
	P.X=x;
	P.Y=y;
	P.Z=GF2X(1);
}

class Kcurve{

	//Koblitz Curve y^2+xy=x^3+a*x^2+1; a=o or 1;

	public:
	mpz_class m;
	GF2X mt; 	
	GF2X a;
	Point G;
	ZZ mu;
	Kcurve(){
		memset(this, 0, sizeof(Kcurve));
	}
	void init(const char* filename){
		
		ifstream fp;
		fp.open(filename, ios::in);
		
		fp>>m;

		string temp;

		fp>>temp;
		hextoGF2X(mt,temp);
		temp.clear();

		fp>>temp;
		hextoGF2X(a,temp);
		temp.clear();
		
		if(IsZero(a)){
			mu=-1;
		}
		else{
			mu=1;
		}
		
		fp>>temp;
		hextoGF2X(this->G.x,temp);
		temp.clear();

		fp>>temp;
		hextoGF2X(this->G.y,temp);
		temp.clear();

		fp.close() ;
	}

	void print(){
		cout<<m<<endl;
		GF2X::HexOutput=1;
		cout<<mt<<endl;
		cout<<a<<endl;
		cout<<G.x<<endl;
		cout<<G.y<<endl;
		GF2X::HexOutput=0;
	}

};

//Preprocessing using solinas parts mod operation

//constant values
int u;//Mu value
complex<mpf_class> tau=(complex<mpf_class>(u,0)+sqrt(complex<mpf_class>(-7,0)))/complex<mpf_class>(2,0);//Tau value


//Round function on mpf_class
mpz_class round(mpf_class val){
	mpf_class r=val-floor(val);
	mpz_class res=static_cast<mpz_class>(floor(val));
	if(r>=0.5) res=res+1;
	return res;
}


//Luca's sequences
vector<mpz_class> U{0,1};
vector<mpz_class> V{2,u};
void init_sequence(mpz_class m){
	int idx;
	for(mpz_class i=2;i<=m;i=i+1){
		idx=mpz_get_si(i.get_mpz_t());
		U.push_back(u*U[idx-1]-2*U[idx-2]);
        V.push_back(u*V[idx-1]-2*V[idx-2]);
	}

}


//Rounding funtion proposed by Solinas
void round_util(mpz_class &q0,mpz_class &q1,mpf_class g0,mpf_class g1){
    mpz_class f0=round(g0);
    mpz_class f1=round(g1);
    mpz_class h0=0,h1=0;
    mpf_class r0=g0-f0;
    mpf_class r1=g1-f1;
    mpf_class r=2*r0+u*r1;
    if(r>=1){
        if(r0-3*u*r1<-1) h1=u;
        else h0=1;
    }
    else{
        if(r0+4*u*r1>=2) h1=u;
    }
    if(r<-1){
        if(r0-3*u*r1>=1) h1=-u;
        else h0=-1;
    }
    else{
        if(r0+4*u*r1<-2) h1=-u;
    }
    q0=f0+h0;
    q1=f1+h1;
    return;
}


//Division function
int M;
void division_util(mpz_class &q0,mpz_class &q1,mpz_class &r0,mpz_class &r1,mpz_class c0,mpz_class c1,mpz_class d0,mpz_class d1){
    mpz_class g0=c0*d0+u*c0*d1+2*c1*d1;
    mpz_class g1=c1*d0-c0*d1;
    mpz_class N=d0*d0+u*d0*d1+2*d1*d1;
    //mpz_class N=pow(2,M)+1-V[M];
    mpf_class l0=static_cast<mpf_class>(g0)/N;
    mpf_class l1=static_cast<mpf_class>(g1)/N;
    round_util(q0,q1,l0,l1);
    r0=c0-d0*q0+2*d1*q1;
    r1=c1-d1*q0-d0*q1-u*d1*q1;
    return;
}

void preprocessing(ZZ &R0,ZZ &R1,mpz_class k,const Kcurve &E){
	//initialising Luca's Sequences and Division function
	init_sequence(E.m);
	M=mpz_get_si((E.m).get_mpz_t());

	pair<mpz_class,mpz_class> tau_m{-2*U[M-1],U[M]};//Tau^m
    mpz_class d0=0,d1=0,rem0=0,rem1=0,q0=0,q1=0,r0=0,r1=0;
    division_util(d0,d1,rem0,rem1,tau_m.first-1,tau_m.second,-1,1);// tau^m -1/tau-1
    division_util(q0,q1,r0,r1,k,0,d0,d1);

    R0=conv<ZZ>(mpz_get_str(nullptr,10,r0.get_mpz_t()));
    R1=conv<ZZ>(mpz_get_str(nullptr,10,r1.get_mpz_t()));
    return;
}

//Scalar Multiplication

void pointadd(LDPoint& R,const LDPoint& P,const Point& Q,const Kcurve& E){
	GF2X A,B,C,D;
	GF2X a=E.a;
	GF2X mt;
	mt=E.mt;
	A=( P.Y + (Q.y*SqrMod(P.Z,mt))%mt)%mt;
	B=(P.X+(Q.x*P.Z)%mt)%mt;
	C=(B*P.Z)%mt;
	R.Z=SqrMod(C,mt);
	D=(Q.x*R.Z)%mt;
	R.X=(SqrMod(A,mt) + (C*(A+SqrMod(B,mt)+(a*C))%mt))%mt;
	GF2X temp=(D+R.X)%mt;
	temp=(temp * ((A*C)+R.Z)%mt)%mt;
	R.Y=(temp+((Q.y+Q.x)*SqrMod(R.Z,mt)))%mt;
}

//T-NAF
vector<int> ZZTNAF(ZZ r0,ZZ r1,const ZZ& mu){	
	vector<int> res;
    ZZ val,temp;
	while(r0!=0 or r1!=0){
		if(abs(r0%2)==1){
			val=(r0+2*u*r1)%4;
            if(val<0) val+=4;
            val=2-val;
            r0=r0-val;
		}
		else val=0;
		res.insert(res.begin(),conv<int>(val));
		temp=r0;
        r0=r1+(mu*r0)/ZZ(2);
        r1=-temp/ZZ(2);
	}
    return res;
}

void scalar_mult(Point& R,const Kcurve& E,ZZ r0,ZZ r1){
	vector<int> res=ZZTNAF(r0,r1,E.mu);
	
	GF2X mt=E.mt;

	Point P;
	P.x=E.G.x;
	P.y=E.G.y;
	Point P1;
	P.negate(P1,mt);

	LDPoint Q;
	LDPoint temp;
	
	if(res[0]==1){
		Q.X=P.x;
		Q.Y=P.y;
		Q.Z=1;
	}

	else if(res[0]==-1){
		Q.X=P.x;
		Q.Y=(P.x+P.y)%mt;
		Q.Z=1;
	}
	
	for(int i=1;i<res.size();i++){
		Q.square(mt);
	
		if(res[i]==1){
			pointadd(temp,Q,P,E);
			Q.set(temp);
		}
		else if(res[i]==-1){
			pointadd(temp,Q,P1,E);
			Q.set(temp);
		}
	}
	Q.toPoint(R,mt);
	return;
}

//T3-NAF representation
vector<ZZ> Rem;
void init_rem(){
    for(int i=1;i<32;i=i+1){
        if(i%8!=0){
            Rem.push_back(conv<ZZ>(i));
            Rem.push_back(conv<ZZ>(-i));
        }
    }
    Rem.push_back(ZZ(0));
    return;
}

vector<int> ZZT3NAF(ZZ r0,ZZ r1,const ZZ& mu){	
	init_rem();
	vector<int> res;
    ZZ val,temp1,temp2,r;
    ZZ c0=r0,c1=r1;
	//c0 and c1 not in set of possible remainders
    while(find(Rem.begin(),Rem.end(),c0)==Rem.end() or find(Rem.begin(),Rem.end(),c1)==Rem.end()){
        ZZ t1=(c0)%ZZ(8);
        ZZ t2=(ZZ(2)*mu*c1)%ZZ(8);
        if(t1<0) t1+=8;
        if(t2<0) t2+=8;
        if(t1!=t2){
            r=(5*c0-2*mu*c1)%ZZ(64);
            if(r<0) r+=64;
            r=32-r;
            for(auto &v:Rem){
                ZZ t=(5*v)%ZZ(64);
                if(t<0) t+=64;
                t=32-t;
                if(t==r){
                    val=v;
                    break;
                }
            }
            c0=c0-val;
        }
        else val=0;
		res.insert(res.begin(),conv<int>(val));
		temp1=c0;
        temp2=c1;
        c0=(-3*u*temp1-2*temp2)/ZZ(8);
        c1=(temp1-2*u*temp2)/ZZ(8);
	}
	res.insert(res.begin(),conv<int>(c0));
	res.insert(res.begin(),conv<int>(c1));
    return res;
}

void scalar_mul(Point& R,const Kcurve& E,ZZ r0,ZZ r1){
	
	vector<int> res=ZZT3NAF(r0,r1,E.mu);

	cout<<"T3-NAF representation : ";
	for(auto &v:res){
		cout<<v<<" ";
	}
	cout<<endl;
	cout<<"length of the representation : "<<res.size()<<endl;
	int cnt=count_if(res.begin(),res.end(),[](int &a){return a!=0;});
	cout<<"Number of non-zeroes : "<<cnt<<endl;

	GF2X mt=E.mt;
	Point P,t;
	P.x=E.G.x;
	P.y=E.G.y;

	Point P1;
	P.negate(P1,mt);

	unordered_map<int,Point> pre_comp;
	for(int i=1;i<32;i++){
		scalar_mult(t,E,conv<ZZ>(i),ZZ(0));
		pre_comp[i]=t;
		t.negate(t,mt);
		pre_comp[-i]=t;
	}

	auto start = high_resolution_clock::now();

	LDPoint Q;
	LDPoint temp;
	
	if(res[0]!=0){
		Q.X=pre_comp[res[0]].x;
		Q.Y=pre_comp[res[0]].y;
		Q.Z=1;
	}
	if(res[1]!=0){
		Q.square(mt);
		pointadd(temp,Q,pre_comp[res[1]],E);
		Q.set(temp);
	}
	for(int i=2;i<res.size();i++){
		Q.square(mt);
		Q.square(mt);
		Q.square(mt);
		if(res[i]!=0){
			pointadd(temp,Q,pre_comp[res[i]],E);
			Q.set(temp);
		}
	}
	auto stop = high_resolution_clock::now();
	auto duration=duration_cast<microseconds>(stop-start);
	cout<<"Time taken to execute : "<<duration.count()<<" ms"<<endl;
	Q.toPoint(R,mt);
	return;
}


int main(){
	Kcurve E;
	E.init("../base_points/K163.txt");
	u=conv<int>(E.mu);
	ZZ r0,r1;
	mpz_class k;
	cin>>k;
	preprocessing(r0,r1,k,E);
	Point res;
	//reduced
	scalar_mul(res,E,r0,r1);
	res.print();
	return 0;
}
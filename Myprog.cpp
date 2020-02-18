#include <NTL/ZZ.h>
#include <NTL/BasicThreadPool.h>
#include <NTL/lzz_pXFactoring.h>
#include "EncryptedArray.h"
#include "Myprog.h"
using namespace std;
NTL_CLIENT
static bool noPrint = true;
int CompareC(EncryptedArray ea,FHESecKey secretKey,Ctxt c1,Ctxt c2)
{
  PlaintextArray p1(ea),p2(ea);
  ea.decrypt(c1, secretKey, p1);
  ea.decrypt(c2, secretKey, p2);
  vector<long> a1,a2;
  a1.resize(ea.size());
  a2.resize(ea.size());
  decode(ea,a1,p1);
  decode(ea,a2,p2);
  if(a1[0]<a2[0])
	return -1;
  else if(a1[0]==a2[0])
	return 0;
  else  return 1;
}
int TestIt(long R, long p, long r, long d, long c, long k, long w, 
               long L, long m, const Vec<long>& gens, const Vec<long>& ords,int a,int b)
{
  char buffer[32];
  vector<long> gens1, ords1;
  convert(gens1, gens);
  convert(ords1, ords);
  FHEcontext context(m, p, r, gens1, ords1);
  buildModChain(context, L, c);
  ZZX G;
  G = context.alMod.getFactorsOverZZ()[0];
  FHESecKey secretKey(context);
  const FHEPubKey& publicKey = secretKey;
  secretKey.GenSecKey(); // A +-1/0 secret key
  addSome1DMatrices(secretKey); // compute key-switching matrices that we need

  EncryptedArray ea(context, G);

  PlaintextArray pA(ea),pB(ea);
  int A=a;
  int B=b;
  encode(ea,pA,long(A));
  encode(ea,pB,long(B));
  Ctxt cA(publicKey), cB(publicKey);
  ea.encrypt(cA, publicKey, pA);
  // {ZZX ppp0; ea.encode(ppp0, p0); c0.DummyEncrypt(ppp0);} // dummy encryption
  ea.encrypt(cB, publicKey, pB); // real encryption

  PlaintextArray pk1(ea),pk2(ea),pk3(ea),pk4(ea),pk5(ea);
  encode(ea,pk1,long(0.8));
  encode(ea,pk2,long(15));
  encode(ea,pk3,long(10));
  encode(ea,pk4,long(150000));
  encode(ea,pk5,long(50000));
  Ctxt ck1(publicKey),ck2(publicKey),ck3(publicKey),ck4(publicKey),ck5(publicKey);
  ea.encrypt(ck1, publicKey, pk1);
  ea.encrypt(ck2, publicKey, pk2);
  ea.encrypt(ck3, publicKey, pk3);
  ea.encrypt(ck4, publicKey, pk4);
  ea.encrypt(ck5, publicKey, pk5);
  
  Ctxt cr1(publicKey),cr2(publicKey),cr3(publicKey);
  for (long i = 0; i < R; i++) {
     cA.multiplyBy(ck1);// A*=0.8  
     cB.multiplyBy(ck2);// B*=15
     if(CompareC(ea,secretKey,cA,cB)>0) cr1=cA;
     else cr1=cB;

     cr1.multiplyBy(ck3);// cr1*=10
     if(CompareC(ea,secretKey,cr1,ck4)<0) cr2=cr1;
     else cr2=ck4;

     if(CompareC(ea,secretKey,cr2,ck5)>0) cr3=cr2;
     else cr3=ck5;
  }

  cA.cleanUp();
  cB.cleanUp();
  cr1.cleanUp();
  cr2.cleanUp();
  cr3.cleanUp();

  PlaintextArray pr(ea);
  ea.decrypt(cr3, secretKey, pr);
  vector<long> ar;
  ar.resize(ea.size());
  decode(ea,ar,pr);
  int rs=ar[0];
  return rs;
}

JNIEXPORT jint JNICALL Java_Myprog_helib(JNIEnv *env, jclass jc,jint a,jint b)
{
  bool dry=false;
  long R=1;
  long p=999983;
  long r=1;
  long d=1;
  long c=2;
  long k=80;
  long L=500;
  long s=0;
  long repeat=1;
  long chosen_m=0;
  Vec<long> mvec;
  Vec<long> gens;
  Vec<long> ords;
  long seed=0;
  long nt=1;
  SetSeed(ZZ(seed));
  SetNumThreads(nt);
  long w = 64; // Hamming weight of secret key
  long m = FindM(k, L, c, p, d, s, chosen_m, !noPrint);
  for (long repeat_cnt = 0; repeat_cnt < repeat; repeat_cnt++) {
    return TestIt(R, p, r, d, c, k, w, L, m, gens, ords, a, b);
  }
}


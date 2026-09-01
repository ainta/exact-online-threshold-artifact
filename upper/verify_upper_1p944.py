#!/usr/bin/env python3
"""Exact-rational verifier for E_star <= 243/125 = 1.944.

Trusted certification decisions use only Python integers and fractions.Fraction.
The verifier checks the rational C^2 profile, its analytic tails, a second
carré-du-champ remainder inequality, and the final continuum barrier integral.
"""
from __future__ import annotations
import json, math, sys, time
from fractions import Fraction as F
from pathlib import Path
from typing import List, Tuple

CERT = Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).with_name('estar_upper_1p944_certificate.json')

class VerificationError(RuntimeError):
    pass

def check(condition, message):
    if not condition:
        raise VerificationError(message)

def frac(s):
    if isinstance(s,int): return F(s,1)
    if isinstance(s,float): raise TypeError('floats forbidden')
    return F(s)

def trim(p: List[F]) -> List[F]:
    q=list(p)
    while len(q)>1 and q[-1]==0:q.pop()
    return q

def add(a,b):
    n=max(len(a),len(b));c=[F(0) for _ in range(n)]
    for i,x in enumerate(a):c[i]+=x
    for i,x in enumerate(b):c[i]+=x
    return trim(c)
def neg(a):return [-x for x in a]
def sub(a,b):return add(a,neg(b))
def scale(a,c):return trim([c*x for x in a])
def mul(a,b):
    c=[F(0) for _ in range(len(a)+len(b)-1)]
    for i,x in enumerate(a):
        for j,y in enumerate(b):c[i+j]+=x*y
    return trim(c)
def deriv(a):return trim([F(i)*a[i] for i in range(1,len(a))] or [F(0)])
def power(a,n):
    r=[F(1)];b=a
    while n:
        if n&1:r=mul(r,b)
        b=mul(b,b);n//=2
    return r

def power_to_bernstein(a: List[F]) -> List[F]:
    a=trim(a);n=len(a)-1
    if n==0:return [a[0]]
    return [sum((a[k]*F(math.comb(i,k),math.comb(n,k)) for k in range(i+1)),F(0)) for i in range(n+1)]
def split_bernstein(b: List[F]) -> Tuple[List[F],List[F]]:
    rows=[list(b)]
    while len(rows[-1])>1:
        old=rows[-1];rows.append([(old[i]+old[i+1])/2 for i in range(len(old)-1)])
    return [row[0] for row in rows],[row[-1] for row in rows][::-1]
def certify_nonnegative_power(p: List[F],max_depth=40):
    stack=[(power_to_bernstein(p),0)];leaves=0;deep=0;minimum=None
    while stack:
        q,d=stack.pop();deep=max(deep,d)
        if min(q)>=0:
            leaves+=1;m=min(q);minimum=m if minimum is None or m<minimum else minimum;continue
        if max(q)<0:return False,{'reason':'negative hull','depth':d,'lo':min(q),'hi':max(q)}
        if d>=max_depth:return False,{'reason':'depth exhausted','depth':d,'lo':min(q),'hi':max(q)}
        l,r=split_bernstein(q);stack.extend([(r,d+1),(l,d+1)])
    return True,{'leaves':leaves,'depth':deep,'min_leaf':minimum}

def quintic(l,r,y0,y1):
    h=r-l;c0=y0[0];c1=h*y0[1];c2=h*h*y0[2]/2
    R0=y1[0]-c0-c1-c2;R1=h*y1[1]-c1-2*c2;R2=h*h*y1[2]-2*c2
    return [c0,c1,c2,10*R0-4*R1+R2/2,-15*R0+7*R1-R2,6*R0-3*R1+R2/2]
def cubic(l,r,y0,y1):
    h=r-l;y00,m0=y0;y11,m1=y1
    return [y00,h*m0,3*(y11-y00)-h*(2*m0+m1),2*(y00-y11)+h*(m0+m1)]
def left_jet(a,d):return [4*a+d/2,-a+d/16,3*d/128]
def right_jet(A):return [4*A+4*A*A,-15*A-31*A*A,52*A+232*A*A]
def bump_basis(j,degree=9):
    p=[F(0),F(0),F(0),F(1),F(-3),F(3),F(-1)]
    p=mul(p,power([F(-1),F(2)],j))
    return p+[F(0)]*(degree+1-len(p))

def segment_constraints(l,r,g,C,Lam,lam):
    h=r-l;z=[l,h]
    gp=scale(deriv(g),1/h);gpp=scale(deriv(gp),1/h);gppp=scale(deriv(gpp),1/h)
    drift=add(scale(z,F(1,2)),g)
    Lg=add(scale(gpp,F(1,2)),mul(drift,gp))
    q=mul(gp,gp);qp=scale(mul(gp,gpp),F(2));qpp=scale(add(mul(gpp,gpp),mul(gp,gppp)),F(2))
    Lq=add(scale(qpp,F(1,2)),mul(drift,qp))
    rem=sub(scale(g,C),q);Lrem=sub(scale(Lg,C),Lq)
    return {
      'g':g,
      'g+z':add(g,z),
      '-(Lg+g)':neg(add(Lg,g)),
      'Lg+Lambda*g':add(Lg,scale(g,Lam)),
      'C*g-gp^2':rem,
      'Lrem+lambda*rem':add(Lrem,scale(rem,lam)),
    }

def build_profile(cert):
    a=frac(cert['left_tail']['a']);d=frac(cert['left_tail']['d']);A=frac(cert['right_tail']['A'])
    knots=[]
    for k in cert['core']['interior_knots']:
        den=k['jet_den'];knots.append((F(k['z_num'],k['z_den']),[F(k['g_num'],den),F(k['gp_num'],den),F(k['gpp_num'],den)]))
    segments=[('left transition',F(-4),F(-3),cubic(F(-4),F(-3),left_jet(a,d)[:2],knots[0][1][:2]))]
    for i in range(len(knots)-1):
        l,y0=knots[i];r,y1=knots[i+1];segments.append((f'interior {i}',l,r,quintic(l,r,y0,y1)))
    l,y0=knots[-1];r=F(4);co=quintic(l,r,y0,right_jet(A))+[F(0)]*4
    for j,e in enumerate(cert['core']['right_bumps']):co=add(co,scale(bump_basis(j),frac(e)))
    segments.append(('right transition',l,r,co))
    z0=[y for z,y in knots if z==0]
    check(len(z0)==1, f'expected exactly one knot at z=0, found {len(z0)}')
    return segments,z0[0],a,d,A

def verify_core(cert):
    C=frac(cert['caps']['C']);Lam=frac(cert['caps']['Lambda']);lam=frac(cert['refinement']['lambda_r']);G=frac(cert['caps']['g0'])
    segments,z0,a,d,A=build_profile(cert)
    gamma,gp0,_=z0
    check(0<gamma<G, f'invalid g(0) bounds: gamma={gamma}, G={G}')
    r0=C*gamma-gp0*gp0;R0=frac(cert['refinement']['r0_lower'])
    check(r0>R0>0, f'invalid remainder bound: r0={r0}, R0={R0}')
    stats={};start=time.time()
    for name,l,r,g in segments:
        for cname,p in segment_constraints(l,r,g,C,Lam,lam).items():
            ok,st=certify_nonnegative_power(p)
            if not ok:raise VerificationError(f'{name}: {cname}: {st}')
            old=stats.get(cname,{'max_depth':0,'leaves':0});old['max_depth']=max(old['max_depth'],st['depth']);old['leaves']+=st['leaves'];stats[cname]=old

    # g is globally C^1 and piecewise C^3.  The generalized Ito formula for
    # rem=Cg-(g')^2 has a local-time term at any jump of rem'.  Verify every
    # such jump is nonnegative, so it can only strengthen Lrem>=-lambda rem.
    def endpoint_jet(seg,side):
        _,l,r,p=seg;h=r-l;x=F(0) if side=='left' else F(1)
        def val(q):return sum(c*x**i for i,c in enumerate(q))
        p1=deriv(p);p2=deriv(p1)
        return val(p),val(p1)/h,val(p2)/(h*h)
    jumps=[]
    # left analytic tail to first polynomial segment
    left=(4*a+d/F(2),-a+d/F(16),3*d/F(128))
    right=endpoint_jet(segments[0],'left')
    check(left[:2]==right[:2], f'C1 mismatch at left interface: {left[:2]} != {right[:2]}')
    jumps.append((-F(4),-2*left[1]*(right[2]-left[2])))
    # all polynomial interfaces
    for i in range(len(segments)-1):
        left=endpoint_jet(segments[i],'right');right=endpoint_jet(segments[i+1],'left')
        check(left[:2]==right[:2], f'C1 mismatch at polynomial interface {i}: {left[:2]} != {right[:2]}')
        jumps.append((segments[i][2],-2*left[1]*(right[2]-left[2])))
    # last polynomial segment to analytic right tail
    left=endpoint_jet(segments[-1],'right');right=tuple(right_jet(A))
    check(left[:2]==right[:2], f'C1 mismatch at right interface: {left[:2]} != {right[:2]}')
    jumps.append((F(4),-2*left[1]*(right[2]-left[2])))
    if any(j<0 for _,j in jumps):raise VerificationError(('negative local-time jump',jumps))
    nonzero=[(x,j) for x,j in jumps if j]
    return {'segments':len(segments),'g0':gamma,'r0':r0,'r0_lower':R0,
            'nonnegative_local_time_jumps':nonzero,'stats':stats,'seconds':time.time()-start}

def verify_left_tail(cert):
    a=frac(cert['left_tail']['a']);d=frac(cert['left_tail']['d']);C=frac(cert['caps']['C']);Lam=frac(cert['caps']['Lambda']);lam=frac(cert['refinement']['lambda_r']);delta=Lam-1
    # y=x^{-1/2}=r/2 with x=-z>=4 and 0<=r<=1.
    y=[F(0),F(1,2)];y2=power(y,2);y3=power(y,3);y6=power(y,6);y7=power(y,7)
    P=add([8*a*a-12*a],add(scale(y3,4*a*d-6*d),add(scale(y6,-4*d*d),scale(y7,-3*d))))
    Q=add(neg(P),scale(add([a],scale(y3,d)),8*delta))
    H=sub(scale(add([a],scale(y3,d)),C),mul(y2,power(add([-a],scale(y3,d/2)),2)))
    # y^2 (Lr+lambda r), written as a polynomial in y and then evaluated at y=r/2.
    coeff=[F(0)]*13
    coeff[0]=C*a*(lam+F(1,2)-a)
    coeff[2]=-a*a*lam
    coeff[3]=C*d*(lam-a/F(2)-F(1,4))
    coeff[5]=a*d*(F(3,2)*a+lam-F(3,4))
    coeff[6]=C*d*d/F(2)
    coeff[7]=F(3,8)*C*d
    coeff[8]=d*d*(F(3,4)*a-lam/F(4)+F(3,8))
    coeff[9]=F(15,8)*a*d
    coeff[11]=-F(3,4)*d**3
    coeff[12]=-F(3,2)*d*d
    E=[c/F(2)**k for k,c in enumerate(coeff)]
    for name,p in [('-(Lg+g)',P),('Lg+Lambda*g',Q),('C*g-gp^2',H),('Lrem+lambda*rem',E)]:
        ok,st=certify_nonnegative_power(p)
        if not ok:raise VerificationError(f'left tail {name}: {st}')
    check(a>1 and d>0, f'invalid left-tail parameters: a={a}, d={d}')
    return {'a':a,'d':d}

def verify_right_tail(cert):
    A=frac(cert['right_tail']['A']);C=frac(cert['caps']['C']);Lam=frac(cert['caps']['Lambda']);lam=frac(cert['refinement']['lambda_r']);delta=Lam-1
    check(A>0, f'invalid right-tail parameter A={A}')
    # Original profile inequalities.
    rho_extra=A/2+46*A*A+31*A**3
    c_bound=64*A+256*A*A+256*A**3
    check(rho_extra<=delta and c_bound<=C,
          f'right-tail cap failed: rho_extra={rho_extra}, delta={delta}, c_bound={c_bound}, C={C}')
    # For E=Lr+lambda*r, with w=A exp(-(z^2-16)/2), exact differentiation gives
    # E=C(lambda-1)zw + sum_{k=2}^6 p_k(z)w^k.
    polys={
      2:{6:F(-1),4:F(9)-lam,2:F(-15)+2*lam,1:C*(lam-F(3,2)),0:F(3)-lam},
      3:{6:F(-10),4:F(56)-4*lam,3:-3*C,2:F(-63)+6*lam,1:2*C,0:F(9)-2*lam},
      4:{6:F(-10),4:F(46)-4*lam,3:-2*C,2:F(-42)+4*lam,1:C,0:F(6)-lam},
      5:{6:F(28),4:F(-66),2:F(30)},
      6:{6:F(16),4:F(-32),2:F(12)},
    }
    tail=F(0)
    for k,p in polys.items():
        for m,c in p.items():
            factor=F(1,4) if m==0 else F(4)**(m-1)
            tail+=abs(c)*factor*A**(k-1)
    lead=C*(lam-1)
    check(tail<lead, f'right-tail remainder failed: tail={tail}, lead={lead}')
    return {'rho_extra_bound':rho_extra,'C_bound':c_bound,'remainder_ratio':tail/lead}

def sqrt_lower(x:F,Q:int)->int:
    check(x>0, f'square-root input must be positive, got {x}')
    return math.isqrt((x.numerator*Q*Q)//x.denominator)

def verify_barrier(cert):
    B=frac(cert['caps']['B']);G=frac(cert['caps']['g0']);C=frac(cert['caps']['C']);Lam=frac(cert['caps']['Lambda'])
    lam=frac(cert['refinement']['lambda_r']);R=frac(cert['refinement']['r0_lower']);qplus=frac(cert['refinement']['positive_power_upper']);qminus=frac(cert['refinement']['negative_power_lower'])
    N=int(cert['refinement']['barrier_subintervals']);bits=int(cert['refinement']['sqrt_bits']);Q=1<<bits
    _,z0,_,_,_=build_profile(cert);gamma=z0[0]
    check(qplus>=2*Lam-1 and qplus==F(25,24),
          f'invalid positive exponent bound: qplus={qplus}, 2Lambda-1={2*Lam-1}')
    check(qminus<=lam-1 and qminus==F(3,2),
          f'invalid negative exponent bound: qminus={qminus}, lambda-1={lam-1}')
    K=R/(lam-2)
    A0=B-C*G
    A1=-(G*G-C*G-K)
    gamma2=gamma*gamma
    check(A0>0 and A1<0 and K>0, f'invalid barrier signs: A0={A0}, A1={A1}, K={K}')
    # Lower envelope after s=x^24:
    # P(x)=A0+A1 x^24+gamma^2 x^25-K x^36.
    # On [a,b], use A0+A1 b^24+gamma^2 a^25-K b^36.
    N24=N**24;N25=N24*N;N36=N**36;total=0;minP=None
    for i in range(N):
        p=A0+A1*F((i+1)**24,N24)+gamma2*F(i**25,N25)-K*F((i+1)**36,N36)
        if not p>0:raise VerificationError(('nonpositive barrier cell',i,p))
        minP=p if minP is None or p<minP else minP
        total+=((i+1)**24-i**24)*sqrt_lower(p,Q)
    lower=F(total,N24*Q);target=F(2,3)*G
    if not lower>target:raise VerificationError(('barrier integral',lower,target))
    return {'lower_integral':lower,'target':target,'margin':lower-target,'min_cell_P':minP,'N':N,'bits':bits}

def main():
    cert=json.loads(CERT.read_text(),parse_float=lambda _:(_ for _ in ()).throw(ValueError('float forbidden')))
    check(cert.get('format')=='estar-rational-profile-v2',
          f"unexpected certificate format: {cert.get('format')!r}")
    print('certificate:',CERT)
    print('core: PASS',verify_core(cert))
    print('left tail: PASS',verify_left_tail(cert))
    print('right tail: PASS',verify_right_tail(cert))
    print('barrier: PASS',verify_barrier(cert))
    print('VERIFIED: E_star <=',cert['caps']['B'])
if __name__=='__main__':main()

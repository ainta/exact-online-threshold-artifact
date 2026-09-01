#!/usr/bin/env python3
from __future__ import annotations
import json,re
from functools import reduce
from operator import mul
from pathlib import Path
ROOT=Path(__file__).resolve().parent

def arr(src,name):
 m=re.search(rf"\b{name}\s*=\s*\{{(.*?)\}}\s*;",src,re.S)
 if not m: raise RuntimeError(f"missing array {name}")
 return [int(x) for x in re.findall(r"(\d+)ULL",m.group(1))]
def const(src,name):
 m=re.search(rf"static constexpr (?:int|u64)\s+{name}\s*=\s*(\d+)(?:ULL)?\s*;",src)
 if not m: raise RuntimeError(f"missing constant {name}")
 return int(m.group(1))
def req(c,m):
 if not c: raise RuntimeError(m)
def main():
 cert=json.loads((ROOT/'certificate.json').read_text()); src=(ROOT/'e_star_lower_bound_verifier.cpp').read_text()
 d=arr(src,'dnum'); w=arr(src,'wnum'); data=cert['weighted_HJB_certificate']; prm=cert['verifier_parameters']; st=cert['statement']['E_star_strict_lower_bound']
 req(d==data['length_numerators'],'dnum mismatch'); req(w==data['unnormalized_weights'],'wnum mismatch')
 C={x:const(src,x) for x in ['PREC','M','D_DEN','GRID_DEN','RADIUS','SCALE_BITS','TARGET_NUM','TARGET_DEN']}
 req(C['M']==data['number_of_intervals']==len(d),'interval count mismatch'); req(C['D_DEN']==data['length_denominator'],'denominator mismatch')
 req(sum(d)==C['D_DEN'],'length sum mismatch'); req(all(x<y for x,y in zip(w,w[1:])),'weights not increasing')
 req(C['TARGET_NUM']==st['numerator'] and C['TARGET_DEN']==st['denominator'],'target mismatch')
 req(C['PREC']==prm['mpfr_precision_bits'],'precision mismatch'); req(C['GRID_DEN']==prm['grid_denominator'],'grid mismatch')
 req(C['RADIUS']==prm['domain_radius'],'radius mismatch'); req(2*C['RADIUS']*C['GRID_DEN']==prm['number_of_grid_intervals'],'interval-grid mismatch')
 req(C['SCALE_BITS']==prm['dyadic_scale_bits'],'scale mismatch')
 block=re.search(r"PRIMES\s*=\s*\{\{\{(.*?)\}\}\}\s*;",src,re.S)
 primes=[int(x) for x in re.findall(r"(\d+)u,\s*\d+u",block.group(1))]
 req(primes==prm['ntt_primes'],'prime mismatch'); crt=reduce(mul,primes,1); req(crt==int(prm['crt_modulus']),'CRT mismatch')
 scale=1<<C['SCALE_BITS']; worst=2*prm['number_of_grid_intervals']*scale*scale
 req(worst==int(prm['worst_case_integer_convolution_bound']),'bound mismatch'); req(worst<crt,'CRT too small')
 req(sum(x*y for x,y in zip(d,w))==int(data['normalization_product_sum']),'normalization mismatch')
 print(f"PASS: certificate.json matches C++ arrays and metadata ({len(d)} stages, {prm['number_of_grid_intervals']} grid intervals)")
if __name__=='__main__': main()

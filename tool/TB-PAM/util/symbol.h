/*!
 \file symbol.h
 \brief 元素記号クラス
*/

#ifndef __SYMBOL_H_INCLUDED
#define __SYMBOL_H_INCLUDED

// 元素名を元素配列のインデックスに変換する構造体
struct ElementSymbol
{
  enum {
    X=0, // unspecified element
     H, He,
    Li, Be,  B,  C,  N,  O,  F, Ne,
    Na, Mg, Al, Si,  P,  S, Cl, Ar,
     K, Ca, Sc, Ti,  V, Cr, Mn, Fe, Co, Ni, Cu, Zn, Ga, Ge, As, Se, Br, Kr,
    Rb, Sr,  Y, Zr, Nb, Mo, Tc, Ru, Rh, Pd, Ag, Cd, In, Sn, Sb, Te,  I, Xe,
    Cs, Ba,
    La, Ce, Pr, Nd, Pm, Sm, Eu, Gd, Tb, Dy, Ho, Er, Tm, Yb, Lu,
    Hf, Ta,  W, Re, Os, Ir, Pt, Au, Hg, Tl, Pb, Bi, Po, At, Rn,
    Fr, Ra,
    Ac, Th, Pa, U, Np, Pu, Am, Cm, Bk, Cf, Es, Fm, Md, No, Lr
  } index;

  inline ElementSymbol( void ){
    index = X;
  }

  ElementSymbol( const char* symbol );

  inline operator unsigned int () const {
    return (unsigned int)index;
  }
};

inline ElementSymbol::ElementSymbol( const char* symbol )
{
  if( symbol[0] == '\0' ){ // in case of null string
    index = X;
    return;
  }
  if( symbol[1] != '\0' && symbol[2] != '\0' ){ // in case of too long string
    index = X;
    return;
  }
  
  switch( symbol[0] ){
  case 'A' : {
    if( false );
    else if( symbol[1] == "Ac"[1] ) index = Ac;
    else if( symbol[1] == "Ag"[1] ) index = Ag;
    else if( symbol[1] == "Al"[1] ) index = Al;
    else if( symbol[1] == "Am"[1] ) index = Am;
    else if( symbol[1] == "Ar"[1] ) index = Ar;
    else if( symbol[1] == "As"[1] ) index = As;
    else if( symbol[1] == "At"[1] ) index = At;
    else if( symbol[1] == "Au"[1] ) index = Au;
    else                            index =  X; // unknown element
  } break;
  case 'B' : {
    if( false );
    else if( symbol[1] == "B"[1]  ) index =  B;
    else if( symbol[1] == "Ba"[1] ) index = Ba;
    else if( symbol[1] == "Be"[1] ) index = Be;
    else if( symbol[1] == "Bi"[1] ) index = Bi;
    else if( symbol[1] == "Bk"[1] ) index = Bk;
    else if( symbol[1] == "Br"[1] ) index = Br;
    else                            index =  X; // unknown element
  } break;
  case 'C' : {
    if( false );
    else if( symbol[1] == "C"[1]  ) index =  C;
    else if( symbol[1] == "Ca"[1] ) index = Ca;
    else if( symbol[1] == "Cd"[1] ) index = Cd;
    else if( symbol[1] == "Ce"[1] ) index = Ce;
    else if( symbol[1] == "Cf"[1] ) index = Cf;
    else if( symbol[1] == "Cl"[1] ) index = Cl;
    else if( symbol[1] == "Cm"[1] ) index = Cm;
    else if( symbol[1] == "Co"[1] ) index = Co;
    else if( symbol[1] == "Cr"[1] ) index = Cr;
    else if( symbol[1] == "Cs"[1] ) index = Cs;
    else if( symbol[1] == "Cu"[1] ) index = Cu;
    else                            index =  X; // unknown element
  } break;
  case 'D' : {
    if( false );
    else if( symbol[1] == "Dy"[1] ) index = Dy;
    else                            index =  X; // unknown element
  } break;
  case 'E' : {
    if( false );
    else if( symbol[1] == "Er"[1] ) index = Er;
    else if( symbol[1] == "Es"[1] ) index = Es;
    else if( symbol[1] == "Eu"[1] ) index = Eu;
    else                            index =  X; // unknown element
  } break;
  case 'F' : {
    if( false );
    else if( symbol[1] == "F"[1]  ) index =  F;
    else if( symbol[1] == "Fe"[1] ) index = Fe;
    else if( symbol[1] == "Fm"[1] ) index = Fm;
    else if( symbol[1] == "Fr"[1] ) index = Fr;
    else                            index =  X; // unknown element
  } break;
  case 'G' : {
    if( false );
    else if( symbol[1] == "Ga"[1] ) index = Ga;
    else if( symbol[1] == "Gd"[1] ) index = Gd;
    else if( symbol[1] == "Ge"[1] ) index = Ge;
    else                            index =  X; // unknown element
  } break;
  case 'H' : {
    if( false );
    else if( symbol[1] == "H"[1]  ) index =  H;
    else if( symbol[1] == "He"[1] ) index = He;
    else if( symbol[1] == "Hf"[1] ) index = Hf;
    else if( symbol[1] == "Hg"[1] ) index = Hg;
    else if( symbol[1] == "Ho"[1] ) index = Ho;
    else                            index =  X; // unknown element
  } break;
  case 'I' : {
    if( false );
    else if( symbol[1] == "I"[1]  ) index =  I;
    else if( symbol[1] == "In"[1] ) index = In;
    else if( symbol[1] == "Ir"[1] ) index = Ir;
    else                            index =  X; // unknown element
  } break;
  case 'J' : {
    if( false );
    else                            index =  X; // unknown element
  } break;
  case 'K' : {
    if( false );
    else if( symbol[1] == "K"[1]  ) index =  K;
    else if( symbol[1] == "Kr"[1] ) index = Kr;
    else                            index =  X; // unknown element
  } break;
  case 'L' : {
    if( false );
    else if( symbol[1] == "La"[1] ) index = La;
    else if( symbol[1] == "Li"[1] ) index = Li;
    else if( symbol[1] == "Lr"[1] ) index = Lr;
    else if( symbol[1] == "Lu"[1] ) index = Lu;
    else                            index =  X; // unknown element
  } break;
  case 'M' : {
    if( false );
    else if( symbol[1] == "Md"[1] ) index = Md;
    else if( symbol[1] == "Mg"[1] ) index = Mg;
    else if( symbol[1] == "Mn"[1] ) index = Mn;
    else if( symbol[1] == "Mo"[1] ) index = Mo;
    else                            index =  X; // unknown element
  } break;
  case 'N' : {
    if( false );
    else if( symbol[1] == "N"[1]  ) index =  N;
    else if( symbol[1] == "Na"[1] ) index = Na;
    else if( symbol[1] == "Nb"[1] ) index = Nb;
    else if( symbol[1] == "Nd"[1] ) index = Nd;
    else if( symbol[1] == "Ne"[1] ) index = Ne;
    else if( symbol[1] == "Ni"[1] ) index = Ni;
    else if( symbol[1] == "No"[1] ) index = No;
    else if( symbol[1] == "Np"[1] ) index = Np;
    else                            index =  X; // unknown element
  } break;
  case 'O' : {
    if( false );
    else if( symbol[1] == "O"[1]  ) index =  O;
    else if( symbol[1] == "Os"[1] ) index = Os;
    else                            index =  X; // unknown element
  } break;
  case 'P' : {
    if( false );
    else if( symbol[1] == "P"[1]  ) index =  P;
    else if( symbol[1] == "Pa"[1] ) index = Pa;
    else if( symbol[1] == "Pb"[1] ) index = Pb;
    else if( symbol[1] == "Pd"[1] ) index = Pd;
    else if( symbol[1] == "Pm"[1] ) index = Pm;
    else if( symbol[1] == "Po"[1] ) index = Po;
    else if( symbol[1] == "Pr"[1] ) index = Pr;
    else if( symbol[1] == "Pt"[1] ) index = Pt;
    else if( symbol[1] == "Pu"[1] ) index = Pu;
    else                            index =  X; // unknown element
  } break;
  case 'Q' : {
    if( false );
    else                            index =  X; // unknown element
  } break;
  case 'R' : {
    if( false );
    else if( symbol[1] == "Ra"[1] ) index = Ra;
    else if( symbol[1] == "Rb"[1] ) index = Rb;
    else if( symbol[1] == "Re"[1] ) index = Re;
    else if( symbol[1] == "Rh"[1] ) index = Rh;
    else if( symbol[1] == "Rn"[1] ) index = Rn;
    else if( symbol[1] == "Ru"[1] ) index = Ru;
    else                            index =  X; // unknown element
  } break;
  case 'S' : {
    if( false );
    else if( symbol[1] == "S"[1]  ) index =  S;
    else if( symbol[1] == "Sb"[1] ) index = Sb;
    else if( symbol[1] == "Sc"[1] ) index = Sc;
    else if( symbol[1] == "Se"[1] ) index = Se;
    else if( symbol[1] == "Si"[1] ) index = Si;
    else if( symbol[1] == "Sm"[1] ) index = Sm;
    else if( symbol[1] == "Sn"[1] ) index = Sn;
    else if( symbol[1] == "Sr"[1] ) index = Sr;
    else                            index =  X; // unknown element
  } break;
  case 'T' : {
    if( false );
    else if( symbol[1] == "Ta"[1] ) index = Ta;
    else if( symbol[1] == "Tb"[1] ) index = Tb;
    else if( symbol[1] == "Tc"[1] ) index = Tc;
    else if( symbol[1] == "Te"[1] ) index = Te;
    else if( symbol[1] == "Th"[1] ) index = Th;
    else if( symbol[1] == "Ti"[1] ) index = Ti;
    else if( symbol[1] == "Tl"[1] ) index = Tl;
    else if( symbol[1] == "Tm"[1] ) index = Tm;
    else                            index =  X; // unknown element
  } break;
  case 'U' : {
    if( false );
    else if( symbol[1] == "U"[1]  ) index =  U;
    else                            index =  X; // unknown element
  } break;
  case 'V' : {
    if( false );
    else if( symbol[1] == "V"[1]  ) index =  V;
    else                            index =  X; // unknown element
  } break;
  case 'W' : {
    if( false );
    else if( symbol[1] == "W"[1]  ) index =  W;
    else                            index =  X; // unknown element
  } break;
  case 'X' : {
    if( false );
    else if( symbol[1] == "X"[1]  ) index =  X; // unknown element
    else if( symbol[1] == "Xe"[1] ) index = Xe;
    else                            index =  X; // unknown element
  } break;
  case 'Y' : {
    if( false );
    else if( symbol[1] == "Y"[1]  ) index =  Y;
    else if( symbol[1] == "Yb"[1] ) index = Yb;
    else                            index =  X; // unknown element
  } break;
  case 'Z' : {
    if( false );
    else if( symbol[1] == "Zn"[1] ) index = Zn;
    else if( symbol[1] == "Zr"[1] ) index = Zr;
    else                            index = X; // unknown element
  } break;
  default : {
    index = X; // unknown element
  } break;
  }
}


#endif // __SYMBOL_H_INCLUDED

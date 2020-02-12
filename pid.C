float FitFunc(float pt,float a0,float a1,float a2,float a3){
return a0+a1*pt+a2/pt+a3*log(pt);
}

static float IsPionW(float m2tof, float pt){
  float ispion=-9999.0;
  float sigma=1.0;
  float m2exp=0.01948816;
  sigma=FitFunc(pt,-0.0269204,0.0222748,0.0180227,0.0229666);
  ispion=(m2tof-m2exp)/sigma;
  return  ispion;
}

static float IsPionE(float m2tof, float pt){
  float ispion=-9999.0;
  float sigma=1.0;
  float m2exp=0.01948816;
  sigma=FitFunc(pt,-0.0430898,0.0620022,0.000438084,-0.0237895);
  ispion=(m2tof-m2exp)/sigma;
  return  ispion;
}

static float IsKaonW(float m2tof, float pt){
  float iskaon=-9999.0;
  float sigma=1.0;
  float m2exp=0.24371698;
  sigma=FitFunc(pt,-0.0752283,0.130887,-0.033046,-0.131315);
  iskaon=(m2tof-m2exp)/sigma;
  return  iskaon;
}

static float IsKaonE(float m2tof, float pt){
  float iskaon=-9999.0;
  float sigma=1.0;
  float m2exp=0.24371698;
  sigma=FitFunc(pt,-0.0963886,0.173101,-0.0522655,-0.18556);
  iskaon=(m2tof-m2exp)/sigma;
  return  iskaon;
}

static float IsProtonW(float m2tof, float pt){
  float iskaon=-9999.0;
  float sigma=1.0;
  float m2exp=0.880;
  sigma=FitFunc(pt,-0.158353,-0.156652,0.371712,0.556339);
  iskaon=(m2tof-m2exp)/sigma;
  return  iskaon;
}

static float IsProtonE(float m2tof, float pt){
  float iskaon=-9999.0;
  float sigma=1.0;
  float m2exp=0.880;
  sigma=FitFunc(pt,-0.177535,-0.0653232,0.300963,0.383232);
  iskaon=(m2tof-m2exp)/sigma;
  return  iskaon;
}

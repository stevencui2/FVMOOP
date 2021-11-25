
int N=7,M=7;double length=1.0;
Grid mygrid_(N,M,length);
int NI,NJ;
NI=mygrid_.pNI();
NJ=mygrid_.pNJ();
Fields fieldsOper(NI,NJ);
Solution sol;
Fields::vectorfields U(NI,Fields::vec1dfields(NJ));
Fields::vectorfields V(NI,Fields::vec1dfields(NJ));
Fields::vectorfields P(NI,Fields::vec1dfields(NJ));
Fields::vectorfields UO(NI,Fields::vec1dfields(NJ));
Fields::vectorfields UOO(NI,Fields::vec1dfields(NJ));
Fields::vectorfields VO(NI,Fields::vec1dfields(NJ));
Fields::vectorfields VOO(NI,Fields::vec1dfields(NJ));
Fields::vectorfields PP(NI,Fields::vec1dfields(NJ));
Fields::vectorfields DPX(NI,Fields::vec1dfields(NJ));
Fields::vectorfields DPY(NI,Fields::vec1dfields(NJ));
Fields::vectorfields massFluxE(NI,Fields::vec1dfields(NJ));
Fields::vectorfields massFluxN(NI,Fields::vec1dfields(NJ));


fieldsOper.getGridInfoPassed(U,mygrid_,sol);
fieldsOper.getGridInfoPassed(V,mygrid_,sol);
fieldsOper.getGridInfoPassed(P,mygrid_,sol);
fieldsOper.getGridInfoPassed(PP,mygrid_,sol);
fieldsOper.getGridInfoPassed(DPX,mygrid_,sol);
fieldsOper.getGridInfoPassed(DPY,mygrid_,sol);
fieldsOper.getGridInfoPassed(massFluxE,mygrid_,sol);
fieldsOper.getGridInfoPassed(massFluxN,mygrid_,sol);

string east="East";
string west="West";
string north="North";
string south="South";

fieldsOper.inletboundaryCondition(U,east,0.0);
fieldsOper.inletboundaryCondition(U,west,0.0);
fieldsOper.inletboundaryCondition(U,north,0.0);
fieldsOper.inletboundaryCondition(U,south,0.0);





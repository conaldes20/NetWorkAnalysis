/*
                        MAIN PROGRAM
        for pipe network analysis by non-linear method using
                        Darcy-Weisbach equation
        Friction factor obtained by iterating the Colebrook-White formula       
*/


#include <iostream>
#include <fstream>
//#include <stdlib>
#include <cstdio>
#include <cmath>
//#include <conio>

const double reo = 997.10;
const double visc = 0.000894;
const double ck = 0.00025908;

void writeToDatFiles(double plRead[],double pdRead[],const int nenw,const int noinw); 
void readData(double HD[], double HDS[], double CSP[], double CSP2[], double CSP3[], int elnu[], int bno[],
        int eno[], double epl[],double epd[], double eff[], const int nenw, const int noinw);
void initialKcoefs(const int cnt, const int elnu[], const int bno[], const int eno[], const double epl[],
        const double epd[], const double H[], double eff[], double kel[],const int nenw);
double fricfac(const double reyn, const double dd);
void KMatrix(const double H[], const double kel[], double ncoef[], double KVals[], const int noinw);
double getKVal(const double kel[], const int elno);
void initialiseNCoef(double ncoef[],const int noinw);
void insertKVal(double ncoef[], const double vok, const int onc);
void insertKsum(double ncoef[], const double sok, const int nonu);
void initialiseKVals(double KVals[],const int noinw);
void initNodalDouble(double PARAM[],const int noinw);
void initNodalInt(int P[],const int noinw);
void initElemDVal(double PA[],const int nenw);
void initElemIVal(int PP[],const int nenw);
void transferncoef(const double ncoef[], double KVals[], const int nonu,const int noinw);
void applyBConditions(double HDS[], double CSP[], double KVals[], int relRC[], int nodeOfKHD[],const int NPH, const int noinw); 
void reducedMatrix(double NKVals[], const double KVals[], const int relRC[], const int UHD,const int noinw);
void reducedHDCSP(double NHDS[], double NCSP[], const double HDS[],const double CSP[], const int relRC[],
        int nodeOfUHD[], const int UHD, const int noinw);
void solveSimulEqns(double coefs[], double x[], double c[], const int UHD);

void initpdRead(double pdRead[],const int nenw);
void insertTOpdRead(double pdRead[], const int elem, const double diam);

void initplRead(double plRead[],const int nenw);
void insertTOplRead(double plRead[], const int elem, const double plen);

void readInt(int& inv); 
void readDouble(double& flov);
void NHIncrement(double H[],double HP[], const double NCSP[], const int nodeOfUHD[],const int UHD, const int noinw);
void ResInConsum(double CSP2[], const double CSP3[],const int noinw);
void ResInHead(double HDS[], const double HD[], const int noinw);
void multi(const double HDS[], double CSP[], const double KVals[],const int noinw);
void Results(const double kel[], double NDS[], double MDS[], const double H[],const int nenw,const int noinw);
void OutPut(const double H[], const double MDS[], const double NDS[],const int nenw, const int noinw);
void AllNodalHeads(double H[], const double HD[],const int noinw);
void CopyOfKVals(double DKVals[], const double KVals[], const int UHD, const int noinw);
void ArraysSet1 (double* &ncoef,double* &CSP, double* &HDS, double* &HD,
         double* &CSP2, double* &CSP3,double* &NDS,double* &H,double* &HP,const int noinw);
void ArraysSet2 (double* &KVals,double* &DKVals, int* &ndnu, int* &relRC, const int noinw);
void ArraysSet3 (double* &pdRead,double* &plRead,double* &MDS, double* &kel, double* &epl,double* &epd, double* &eff,
         int* &elnu, int* &bno, int* &eno, const int nenw);
void ArraysSet4 (double* &NHDS, double* &NCSP, double* &NKVals, int* &nodeOfKHD,
        int* &nodeOfUHD, const int UHD, const int NPH);
void NetWorkAttributes (int& noinw, int& nenw);
void ReadIteTol(double& tol, int& mnite);
void ReadNWAttributes (int& noinw, int& nenw);
void outputData(const int ndnu[], const double HDS[], const double CSP[], const int elnu[], const int bno[],
        const int eno[], const double epl[], const double epd[], const double eff[], const int nenw, const int noinw);
void writeToFiles78(int& NPH,int& NCS, int& UHD,const int noinw);
void readCHD(double& flov);

void readNode(int& inv, const int noinw);
void readElem(int& inv, const int nenw);

void PrintHeader();
void PrintIteTol(const double tol, const int mnite);
void PrintErrConv(const double error, const int istep,const double tol);
void PrintData(const int ndnu[], const double HDS[], const double CSP[], const int elnu[], const int bno[],
        const int eno[], const double epl[], const double epd[], const double eff[], const int nenw, const int noinw);
void PrintResult(const double H[], const double MDS[], const double NDS[],const int nenw, const int noinw);
  
main() {        
        int mchoice = 0;
        int schoice = 0;
        int nparams = 0;
        int flag = 0;
        int nenw;       //Number of elements in a network
        int noinw;      //Number of nodes in a network
        int NPH;        //Number of prescribed heads
        int NCS;        //Number of specified consumption
        int UHD;        //Number of unknown heads
        int mnite;      //Number of iterations
        double tol;     //Tolerance 
                
        NetWorkAttributes (noinw, nenw);  
        int choice = 0;
                                        
        do {
                std::cout << std::endl << "      PRESCRIBED DATA" << std::endl;
                std::cout << "      1 Enter New Test Data" << std::endl;
                std::cout << "      2 Use Existing Test Data" << std::endl;
                std::cout << "      Enter your choice: ";
                readInt(choice);
                switch(choice) {
                        case 1:
                                writeToFiles78(NPH,NCS, UHD, noinw);
                                break;
                        case 2:
                                std::ifstream inf14("zon114.dat", std::ios::in);
                                while(inf14 >> NPH >> NCS >> UHD) {

                                }
                                inf14.close();
                                break;
                }
        } while((choice != 2) && (choice != 1));
        double* ncoef;
        double* CSP;  //equivalent of al[j]
        double* HDS;  //equivalent of h[j]
        double* HD;   //equivalent of h[j]
        double* CSP2; //equivalent of q[j]
        double* CSP3; //equivalent of elre[j]
        double* NDS;
        double* H;
        double* HP;    //All of size noinw        
        double* KVals;
        double* DKVals; //All of size noinw*noinw
        int* ndnu;
        int* relRC;  //All of size noinw
        double* pdRead;
        double* plRead;
        double* MDS;
        double* kel;
        double* epl;
        double* epd;
        double* eff;  //All of size nenw
        int* elnu;
        int* bno;
        int* eno;    //All of size nenw
        double* NCSP;
        double* NHDS;
        double* NKVals; //All of size UHD*UHD
	double* errtol;
        int* nodeOfKHD;   //Size NPH     Stop here
        int* nodeOfUHD;   //Size UHD
        
        ArraysSet1 (ncoef,CSP, HDS, HD, CSP2, CSP3,NDS,H,HP,noinw);
        ArraysSet2 (KVals,DKVals, ndnu, relRC, noinw);
        ArraysSet3 (pdRead,plRead,MDS, kel, epl,epd,eff,elnu, bno,eno, nenw);
        ArraysSet4 (NHDS,NCSP, NKVals, nodeOfKHD,nodeOfUHD, UHD, NPH);
        //ArraysSet5 (errtol, 2);   
    
        do {
                std::cout << std::endl << "      MAIN MENU" << std::endl;
                std::cout << "      1 Enter New Data for Network" << std::endl;
                std::cout << "      2 Generate Results and Display Output" << std::endl;
                std::cout << "      3 Quit Program" << std::endl;
                std::cout << "      Enter your choice: ";
                readInt(mchoice);                
                if (mchoice == 1) {                               
                                do {
                                        std::cout << std::endl << "      Sub Menu" << std::endl;
                                        std::cout << "      1 Read Network Data" << std::endl;
                                        std::cout << "      2 Back to Main Menu" << std::endl;
                                        std::cout << "      Enter your choice: ";
                                        readInt(schoice);
                                        
                                        switch(schoice) {
                                                case 1:
                                                        writeToDatFiles(plRead,pdRead, nenw, noinw );
                                                        std::ofstream outf6("zon16.dat", std::ios::out);
                                                        nparams = 1;
                                                        outf6 << nparams << std::endl;
                                                        outf6.close();
                                                        break;
                                        }
                                } while(schoice != 2);                                
                } else if(mchoice == 2) {
                                std::ifstream inf6("zon16.dat", std::ios::in);  
                                if(!inf6) {
                                        std::cout << std::endl << "      File could not open" << std::endl;
                                        std::cout << "      Try again!" << std::endl;
                                } else {                                 
                                        inf6 >> nparams;
                                        inf6.close();
                                }
                                
                                std::ifstream inf7("zon17.dat", std::ios::in);
                                std::ifstream inf3("zon13.dat", std::ios::in);
                                std::ifstream inf8("zon18.dat", std::ios::in);
                                if ((!inf7) || (!inf8) || (!inf3)) {
                                    std::cout << std::endl << "     Data Files not existing!" << std::endl;
                                    flag = 0;
                                    break;
                                } else {                                    
                                    flag = 1;
                                    inf7.close();
                                    inf8.close();
                                    inf3.close();
                                }
                                     
                                if((nparams == 1) && (flag == 1)){ 
                                        int count = 0;
                                        int istep = 0;
                                        int again = 0;
                                        ReadIteTol(tol, mnite);
                                        //insert function
                                        initNodalDouble(H,noinw);
                                        initNodalDouble(HP,noinw);
                                        initNodalDouble(HDS,noinw);
                                        initNodalDouble(HD,noinw);
                                        initNodalDouble(CSP,noinw);
                                        initNodalDouble(CSP2,noinw);
                                        initNodalDouble(CSP3,noinw);
                                        //initNodalDouble(x,noinw);
                                        //initNodalDouble(y,noinw);
                                        initNodalDouble(ncoef,noinw);                                   
                                        initNodalInt(ndnu,noinw);
                                        initNodalDouble(NDS,noinw);
                                        initElemDVal(MDS,nenw);
                                        initElemDVal(epd,nenw);
                                        initElemDVal(eff,nenw);
                                        initElemDVal(kel,nenw);
                                        initElemIVal(elnu,nenw);
                                        initElemIVal(bno,nenw);
                                        initElemIVal(eno,nenw);
                                        initialiseKVals(KVals, noinw);
                                        initialiseKVals(DKVals, noinw);                                 
                                        readData(HD, HDS, CSP, CSP2, CSP3, elnu, bno,eno,epl, epd, eff, nenw, noinw);
                                      
                                        //initialiseHP(HP,HDS,noinw);

                                        initialKcoefs(count, elnu, bno,eno, epl,epd, H, eff, kel,nenw);
                                        KMatrix(H, kel, ncoef, KVals,noinw);
                                        CopyOfKVals(DKVals, KVals, UHD, noinw);
                                        applyBConditions(HDS, CSP, KVals, relRC, nodeOfKHD,NPH, noinw);
                                        reducedMatrix(NKVals, KVals, relRC, UHD, noinw);
                                        reducedHDCSP(NHDS, NCSP, HDS,CSP, relRC,nodeOfUHD, UHD, noinw);
                                        solveSimulEqns(NKVals, NHDS, NCSP, UHD);
                                        do {
                                                initialiseKVals(KVals,noinw);
                                                NHIncrement(H, HP,NCSP, nodeOfUHD,UHD, noinw);
                                                count = count + 1;
                                                initialKcoefs(count, elnu, bno,eno, epl,epd, H, eff, kel,nenw);
                                                KMatrix(H, kel, ncoef, KVals,noinw);
                                                CopyOfKVals(DKVals, KVals, UHD, noinw);
                                                ResInConsum(CSP2, CSP3,noinw);
                                                ResInHead(HDS, HD, noinw);
                                                applyBConditions(HDS, CSP, KVals, relRC, nodeOfKHD,NPH, noinw);
                                                multi(H, CSP, KVals,noinw);

                                                double diff = 0.00;
                                                double del = 0.0;
                                                double tot = 0.0;
                                                //double tol = 0.001;
                                                double error = 0.0;
                                                for(int z = 0; z < noinw; z++) {
                                                        CSP[z] = CSP2[z] - CSP[z];
                                                        //del = del + pow(CSP[z], 2);
                                                        //tot = tot + fabs(CSP2[z]);
                                                        diff = HP[z] - H[z];
                                                        del = del + pow(diff, 2);
                                                        tot = tot + fabs(H[z]);
                                                }
                                                if (tot == 0.00) {
                                                        tot = 1.00;                                             }
                                                error = sqrt(del)/tot;
                                                //std::cout << "       " << std::endl;
                                                //std::cout << "      Error = " << error <<  "  Tolerance = " << tol << std::endl;
                                                //int xc;
                                                //std::cout << "     Enter 1";
                                                //std::cin >> xc;
                                                if ((error - tol) <= 0.0) {
                                                        //clrscr();
                                                        std::cout << std::endl << "       " << std::endl;
                                                        std::cout << std::endl << "      Error = " << error <<  "  Tolerance = " << tol << std::endl;
                                                        std::cout << std::endl << "      Convergence achieved on cycle " << istep << std::endl; 
														
                                                        //call result
                                                        AllNodalHeads(H, HD,noinw);                                                     
                                                        outputData(ndnu, HD, CSP3, elnu, bno,eno, epl,epd, eff, nenw, noinw);                                                 
                                                        Results(kel, NDS, MDS, H,nenw,noinw);   
                                                        OutPut(H, MDS, NDS,nenw, noinw);
                                                        again = 0;
                                                        std::cout << std::endl;
                                                        std::cout << std::endl << "     Is a printer ready (Y/N) ?";
                                                        char ans;               
                                                        std::cin >> ans;                                                                                          
                                                        if ((ans == 'Y') || (ans == 'y')) { 
                                                            PrintHeader();
                                                            PrintIteTol(tol, mnite);
                                                            PrintErrConv(error, istep, tol);
                                                            PrintData(ndnu, HDS, CSP, elnu, bno, eno, epl, epd, eff, nenw, noinw);
                                                            PrintResult(H, MDS, NDS,nenw, noinw);
                                                        }  
                                               } else if ((error - tol) > 0.0) {
                                                        istep = istep + 1;
                                                        if ((istep - mnite) > 0) {
                                                                std::cout << "      Maximum number of iteration reached without convergence" << std::endl;
                                                                //call result
                                                                AllNodalHeads(H, HD,noinw);
                                                                outputData(ndnu, HD, CSP3, elnu, bno,eno, epl,epd, eff, nenw, noinw);
                                                                Results(kel, NDS, MDS, H,nenw,noinw);   
                                                                OutPut(H, MDS, NDS,nenw, noinw);
                                                                again = 0;
                                                                std::cout << std::endl;
                                                                std::cout << std::endl << "     Is a printer ready (Y/N) ?";
                                                                char ans;               
                                                                std::cin >> ans;                                                                                          
                                                                if ((ans == 'Y') || (ans == 'y')) { 
                                                                    PrintHeader();
                                                                    PrintIteTol(tol, mnite);
                                                                    PrintErrConv(error, istep, tol);
                                                                    PrintData(ndnu, HDS, CSP, elnu, bno, eno, epl, epd, eff, nenw, noinw);
                                                                    PrintResult(H, MDS, NDS,nenw, noinw);
                                                                }  
                                                        } else if ((istep - mnite) <= 0) {
                                                                //Compute tangent system matrix
                                                                for(int w = 0; w < noinw; w++)
                                                                        CSP[w] = CSP[w]/0.50;
                                                                //Obtain new nodal head increment
                                                                reducedMatrix(NKVals, KVals, relRC, UHD,noinw);
                                                                reducedHDCSP(NHDS, NCSP, HDS,CSP, relRC,nodeOfUHD, UHD, noinw);
                                                                solveSimulEqns(NKVals, NHDS, NCSP, UHD);
                                                                //Back and compute current total heads
                                                                again = 1;
                                                        }
                                                }                                       
                                        } while(again > 0);
                                } else if (nparams == 0) {
                                        std::cout << std::endl << "      Network Data for Computation"
                                        " Not Available" << std::endl;
                                        std::cout << std::endl << "      Back to Main Menu and Pick Option 1" << std::endl;
                                }                                            
                }

        } while(mchoice != 3);
        delete[] ncoef;
        delete[] CSP;
        delete[] HDS;
        delete[] HD;
        delete[] CSP2;
        delete[] CSP3;
        delete[] NDS;
        delete[] H;
        delete[] HP;
        //delete[] x;
        //delete[] y;
        delete[] KVals;
        delete[] DKVals;
        delete[] ndnu;
        delete[] relRC;
        delete[] pdRead;
        delete[] plRead;
        delete[] MDS;
        delete[] kel;
        delete[] epl;
        delete[] epd;
        delete[] eff;
        delete[] elnu;
        delete[] bno;
        delete[] eno;
        delete[] NHDS;        
        delete[] NCSP;
        delete[] NKVals;
        delete[] nodeOfKHD;
        delete[] nodeOfUHD;
	delete[] errtol;
        return 0;
}

void writeToDatFiles(double plRead[],double pdRead[],const int nenw,const int noinw ) {
        int thro = 0;
        int fthro2 = 0;
        int fthro3 = 0;
        int nodnu = 0;
        int nuelem = 0;
        int bnode = 0;
        int enode = 0;
        int elem = 0;
        double pdiam = 0.0;
        double plen = 0.0;
        double ffact = 0.0;
        //double xcoord = 0.0;
        //double ycoord = 0.0;
        std::ofstream outf2("zon12.dat", std::ios::app);
        std::ofstream outf3("zon13.dat", std::ios::app);
        //std::ofstream outf4("zon14.dat", std::ios::out);
        initpdRead(pdRead, nenw);
        initplRead(plRead, nenw);   

        do {
                if (fthro2 == 0) {
                        std::cout << "     Enter a node's number:";
                        fthro2 = 1;
                } else if (fthro2 == 1)
                        std::cout << "     Enter next node's number:";
                readNode(nodnu, noinw);

                //std::cout << "     Enter x-coordinate of this node:";
                //readCHD(xcoord);
                //std::cout << "     Enter y-coordinate of this node:";
                //readCHD(ycoord);

                std::cout << "     How many elements meet at this node?";
                readInt(nuelem);

                outf2 << nodnu << " " << nuelem << std::endl;
                //outf4 << nodnu <<  " " << xcoord << " " << ycoord << std::endl;
                int noel = 0;
                fthro3 = 0;
                bnode = nodnu;
                for(noel = 1; noel <= nuelem; noel++) {
                        if (fthro3 == 0) {
                                std::cout << "     Enter the number of one such elements:";
                                fthro3 = 1;
                        } else if (fthro3 == 1)
                                std::cout << "     Enter the number of another such elements:";
                        readElem(elem, nenw);
                        std::cout << "     Enter a node number this element connects:";
                        readNode(enode, noinw);

                        plen = plRead[(elem - 1)];
                        if (plen == 0.0) {
                                std::cout << "     Enter pipe's length for this element:";
                                readDouble(plen);
                                insertTOplRead(plRead, elem, plen);

                        }
                        
                        pdiam = pdRead[(elem - 1)];
                        if (pdiam == 0.0) {
                                std::cout << "     Enter pipe's diameter for this element:";
                                readDouble(pdiam);
                                insertTOpdRead(pdRead, elem, pdiam);

                        }
                        ffact = 0.005;
                        outf3 << elem <<  " " << bnode << " " << enode << " " << plen << " " << pdiam << " " << ffact  << std::endl;
                        
                }

                char ans;

                std::cout << std::endl << "     Do you want to continue with the data entry (Y/N)? ";
                //std::cout << std::endl;  
                std::cin >> ans;
                                                     
                if ((ans == 'Y') || (ans == 'y')) 
                        thro = 0;
                else {
                        thro = 1;
                        break;
                }            

        } while (thro == 0);
        
        outf2.close();
        outf3.close();
        //outf4.close();  
        std::cout << "   " << std::endl;
}

void writeToFiles78(int& NPH,int& NCS,int& UHD,const int noinw) {
        int cnod = 0;
        int nodnu = 0;  
        double head = 0.0;
        double consum = 0.0;
        int fthro2 = 0;
        std::cout << "   " << std::endl;
        std::cout << "     BOUNDARY DATA: PRESCRIBED HEADS";
        std::cout << std::endl;
        std::cout << "     Enter number of nodes with prescribed heads: ";                          
        readInt(NPH);           
        std::cout << std::endl;
        
        std::ofstream outf7("zon17.dat", std::ios::out);
        for(cnod = 1; cnod <= NPH; cnod++) {
                if (fthro2 == 0) {
                        std::cout << "     Node's number with a prescribed head:";
                        fthro2 = 1;
                } else if (fthro2 == 1)
                        std::cout << "     Next Node's number with prescribed head:";
                readNode(nodnu, noinw);
                std::cout << "     Enter head's value for this node:";
                readCHD(head);
                outf7 << nodnu << " " << head << std::endl;
        }
        outf7.close();
        UHD = noinw - NPH;

        fthro2 = 0;
        std::cout << "   " << std::endl;
        std::cout << "     BOUNDARY DATA: NODAL DISCHARGES";
        std::cout << "   " << std::endl;
        std::cout << "     Enter number of nodes with specified discharges: ";              
        readInt(NCS);
        std::ofstream outf8("zon18.dat", std::ios::out);
        for(cnod = 1; cnod <= NCS; cnod++) {   
                if (fthro2 == 0) {
                        std::cout << "     Node's number with known nodal discharge:";
                        fthro2 = 1;
                } else if (fthro2 == 1)
                        std::cout << "     Next node's number with known nodal discharge:";
                readNode(nodnu, noinw);
                std::cout << "     Enter nodal discarge for this node:";            
                readCHD(consum);
                outf8 << nodnu << " " << consum << std::endl;
        }
        outf8.close();
        std::ofstream outf14("zon114.dat", std::ios::out);
        outf14 << NPH << " " << NCS << " " << UHD << std::endl;
        outf14.close();
}

void readData(double HD[], double HDS[], double CSP[], double CSP2[], double CSP3[], int elnu[], int bno[],
int eno[], double epl[], double epd[], double eff[], const int nenw, const int noinw) {
        int nodnu = 0;
        double head = 0.0;
        double consum = 0.0;
        //double xcoord = 0.0;
        //double ycoord = 0.0;
        int elem = 0;
        int bnode = 0;
        int enode = 0;
        double plen = 0.0;
        double pdiam = 0.0;
        double ffact = 0.0;
        //int x1;
        int n = 0;
        
        std::ifstream inf7("zon17.dat", std::ios::in);
        while(inf7 >> nodnu >> head) {
                //std::cout << "     nodnu:  " << nodnu <<  "    head:  " << head << std::endl;
                for(n = 0; n < noinw; n++) {
                        if ((nodnu - 1) == n) {
                                HDS[n] = head;
                                HD[n] = head;
                                break;
                        }
                }
        }
        inf7.close();
        
        std::ifstream inf8("zon18.dat", std::ios::in);
        while(inf8 >> nodnu >> consum) {
                //std::cout << "     nodnu:  " << nodnu <<  "    consum:  " << consum << std::endl;
                for(n = 0; n < noinw; n++) {
                        if ((nodnu - 1) == n) {
                                CSP[n] = consum;
                                CSP2[n] = consum;
                                CSP3[n] = consum;
                                break;
                        }
                }
        }
        inf8.close();

        n = 0;
        std::ifstream inf3("zon13.dat", std::ios::in);
        while(inf3 >> elem >> bnode >> enode >>  plen >> pdiam >> ffact) {
                //std::cout << "     elem:  " << elem << "  bnode:  " << bnode << "  enode:  " << enode << std::endl;
                for(n = 0; n < nenw; n++) {
                        if ((elem - 1) == n) {
                                elnu[n] = elem;
                                bno[n] = bnode;
                                eno[n] = enode;
                                epl[n] = plen;
                                epd[n] = pdiam;
                                eff[n] = ffact;
                                //std::cout << "     elem:  " << elnu[n] << "  bnode:  " << bno[n] << "  enode:  " << eno[n] <<  "  pdiam:  " << epd[n] << "  pffact:  " << eff[n] << std::endl;
                                break;
                        }
                }
        }
        inf3.close();
}

void initialKcoefs(const int cnt, const int elnu[], const int bno[], const int eno[], const double epl[],
        const double epd[], const double H[], double eff[], double kel[],const int nenw) {
        double plen = 0.0;
        double elr = 0.0;
        double reyn = 0.0;
        double dd = 0.0;
        double e = 0.0;
        int bn = 0;
        int en = 0;
       
        e = 4.0*reo/(3.1416*visc);
        if (cnt == 0) {         
                int n = 0;
                int el = 0;
                for(n = 0; n < nenw; n++) {
                        el = elnu[n];
                        bn = bno[el - 1];
                        en = eno[el - 1];
                        bn = bn - 1;
                        en = en - 1;    
                        //std::cout << "elnu[" << n << "]     bno[" << el - 1 << "]     eno[" << el - 1 << "]  :  " << elnu[n] << "  " << bno[el - 1] << "  " << eno[el - 1] << std::endl;
                        //std::cout << "bno[" << el - 1 << "]     x[" << bn << "]     y[" << bn << "]  :  " << bno[el - 1] << "  " << x[bn] << "  " << y[bn] << std::endl;
                        //std::cout << "eno[" << el - 1 << "]     x[" << en << "]     y[" << en << "]  :  " << eno[el - 1] << "  " << x[en] << "  " << y[en] << std::endl;
                        //std::cout << "pow(epd[n],2.5)   sqrt(plen*eff[n])  :  " << pow(epd[n],2.5) << "  " << sqrt(plen*eff[n]) << std::endl;
                        //std::cout << "plen    kel[" << n << "]  :  " << plen << "  " << kel[n] << std::endl;
                        //std::cout << "     Enter 1";
                        //std::cin >> x1;
                                   
                        plen = epl[n];        //sqrt(pow((x[en] - x[bn]),2.0) + pow((y[en] - y[bn]),2.0));                       
                        kel[n] = 1.7409*(pow(epd[n],2.5))/(sqrt(plen*eff[n]));                  
                }
                
        } else if (cnt > 0) {
                int n = 0;
                int el = 0;
                int a = 0;
                double sgn = 0.0;
                for(n = 0; n < nenw; n++) {
                        el = elnu[n];
                        bn = bno[el];
                        en = eno[el];
                        bn = bno[el - 1];
                        en = eno[el - 1];
                        bn = bn - 1;
                        en = en - 1;                    
                        plen = epl[n];        //sqrt(pow((x[en] - x[bn]),2.0) + pow((y[en] - y[bn]),2.0));
                        //Blasius Equation to compute friction factor
                        a = H[(bn)] - H[(en)];
                        if (a >= 0.0)
                                sgn = 1.0;
                        else
                                sgn = -1.0;
                        elr = sgn*kel[n]*sqrt(fabs(a));
                        reyn = (e*(fabs(elr)))/epd[n];
                        dd = epd[n];
                        if (reyn < 4000.0)
                                reyn = 4000.0;
                        eff[n] = fricfac(reyn, dd);                     
                        if (eff[n] <= 0.0)
                                continue;                       
                        kel[n] = 1.7409*(pow(epd[n],2.5))/(sqrt(plen*eff[n]));
                        //std::cout << "epd[" << n << "]     eff[" << n << "]  :  " << epd[n] << "  " << eff[n] << std::endl;
                        //std::cout << "pow(epd[n],2.5)   sqrt(plen*eff[n])  :  " << pow(epd[n],2.5) << "  " << sqrt(plen*eff[n]) << std::endl;
                        //std::cout << "plen    kel[" << n << "]  :  " << plen << "  " << kel[n] << std::endl;
                       
                }
        }
}

double fricfac(const double reyn, const double dd) {
        double fff = 0.0;
        double ckd = 0.0;
        double vre = 0.0;
        double am = 0.0;
        double cc = 0.0;
        double b = 0.0;
        float jf = 0.0f;
        double jd = 0.0;
        ckd = ck/(3.71*dd);
        vre = 1.26/reyn;
        fff = 0.0001;
        for(int j = 0; j < 400; j++) {
                am = -4.0*(log10(ckd + (vre/(sqrt(fff)))));
                b = 1.0/(am*am);
                cc = (fabs(b - fff))/b;
                jf = float(j);
                jd = double(jf);
                fff = fff + (0.000001*jd);
                if(cc <= 0.01)
                        break;
        }
        return fff;
}

void KMatrix(const double H[], const double kel[], double ncoef[], double KVals[], const int noinw) {
        int nodnu = 0;
        int nuelem = 0;
        int bnode = 0;
        int enode = 0;
        int elem = 0;
        double plen = 0.0;
        double pdiam = 0.0;
        double ffact = 0.0;
        std::ifstream inf2("zon12.dat", std::ios::in);
        std::ifstream inf3("zon13.dat", std::ios::in);

        while(inf2 >> nodnu >> nuelem) {
                int cnt = 0;
                initialiseNCoef(ncoef,noinw);
                //std::cout << "     nodnu:  " << nodnu << "  nuelem:  " << nuelem << std::endl;
                double ksum = 0.0;
                double dh = 0.0;
                double dev = 0.0;
                double kva = 0.0;               
                int nk = 0.0;
                while(inf3 >> elem >> bnode >> enode >>  plen >> pdiam >> ffact) {
                        nk = elem;
                        cnt = cnt + 1;
                        //std::cout << "     elem:  " << elem << "  bnode:  " << bnode << "  enode: " << enode << "  pdiam:  " << pdiam << "  ffact:  " << ffact << std::endl;
                        dh = H[(bnode - 1)] - H[(enode - 1)];                   
                        dev = fabs(dh) - 0.0000001;
                        kva = getKVal(kel, nk);                 
                        if(dev > 0.0)
                                kva = kva*(pow(fabs(dh),-0.50));                        
                        //insert function
                        ksum = ksum + kva;
                        //insert function
                        insertKVal(ncoef, kva, enode);                  
                        if (cnt == nuelem)
                                break;
                }
                //insert function
                insertKsum(ncoef, ksum, nodnu);         
                transferncoef(ncoef, KVals, nodnu,noinw);
        }
        
        inf2.close();
        inf3.close();
        
        //for(int d = 0; d < noinw*noinw; d++) {
        //              std::cout << "KVals[" << d << "]  =  " << KVals[d] << "  ";
        //              std::cout << "     Enter 1";
        //              std::cin >> x1;
        //}       

}


double getKVal(const double kel[], const int elno) {
        double kv = 0.0;
        int n = 0;
        n = elno - 1;
        kv = kel[n];
        return kv;
}


void initialiseNCoef(double ncoef[],const int noinw) {
        for(int n = 0; n < noinw; n++)
                ncoef[n] = 0.0;
}

void insertKVal(double ncoef[], const double vok, const int onc) {
        int n;
        n = onc - 1;    
        ncoef[n] = -vok;
}

void insertKsum(double ncoef[], const double sok, const int nonu) {
        int n;
        n = nonu - 1;   
        ncoef[n] = sok;                 
}

void initialiseKVals(double KVals[],const int noinw) {
        int ind;
        for(int n = 0; n < noinw; n++) {
                for(int m = 0; m < noinw; m++)
                        ind = noinw*n + m;
                        KVals[ind] = 0.0;
        }
}

void initNodalDouble(double PARAM[],const int noinw) {
        for(int m = 0; m < noinw; m++)
                        PARAM[m] = 0.0;
}

void initNodalInt(int P[],const int noinw) {
        for(int m = 0; m < noinw; m++)
                        P[m] = 0;
}       

void initElemDVal(double PA[],const int nenw) {
        for(int m = 0; m < nenw; m++)
                        PA[m] = 0.0;
}

void initElemIVal(int PP[],const int nenw) {
        for(int m = 0; m < nenw; m++)
                        PP[m] = 0;
}

void transferncoef(const double ncoef[], double KVals[], const int nonu,const int noinw) {
        int n;
        int ind;
        
        n = nonu - 1;
        for(int m = 0; m < noinw; m++) {
                        ind = noinw*n + m;
                        KVals[ind] = ncoef[m];
        }
}

void applyBConditions(double HDS[], double CSP[], double KVals[], int relRC[], int nodeOfKHD[],const int NPH, const int noinw) {
        int w = 0;
        int ind;
        for(int t = 0; t < NPH; t++)
                nodeOfKHD[t] = 0;
        for(int h = 0; h < noinw; h++)
                relRC[h] = 1;
        for(int i = 0; i < noinw; i++) {
                if(HDS[i] != 0.0) {
                        nodeOfKHD[w] = i + 1;
                        relRC[i] = 0;
                        for(int j = 0; j < noinw; j++) {
                                ind = noinw*j + i;
                                CSP[j] = CSP[j] - (HDS[i]*KVals[ind]);
                        }
                        for(int jl = 0; jl < noinw; jl++) {
                                if(i != jl) {
                                        ind = noinw*jl + i;
                                        KVals[ind] = 0.0;
                                } else if(i == jl) {
                                        ind = noinw*jl + i;
                                        KVals[ind] = 1.0;
                                }
                        }
                        for(int k = 0; k < noinw; k++) {
                                if(i != k) {
                                        ind = noinw*i + k;
                                        KVals[ind] = 0.0;
                                } else if(i == k) {
                                        ind = noinw*i + k;
                                        KVals[ind] = 1.0;
                                }
                        }
                        w = w + 1;
                }
        }
}

void reducedMatrix(double NKVals[], const double KVals[], const int relRC[], const int UHD,const int noinw) {
        int k = 0;
        int ind, id;
        for(int r = 0; r < UHD; r++) {
                for(int s = 0; s < UHD; s++) {
                        ind = UHD*r + s;
                        NKVals[ind] = 0.0;
                }
        }
        
        for(int n = 0; n < noinw; n++) {
                if(relRC[n] == 1) {
                        int kk = 0;
                        for(int i = 0; i < noinw; i++) {
                                if(relRC[i] == 1) {
                                        ind = UHD*k + kk;
                                        id = noinw*n + i;
                                        NKVals[ind] = KVals[id];
                                        kk = kk + 1;
                                }
                        }
                        k = k + 1;
                }
        }
}

void reducedHDCSP(double NHDS[], double NCSP[], const double HDS[],const double CSP[], const int relRC[],
        int nodeOfUHD[], const int UHD, const int noinw) {      
        int k = 0;
        for(int r = 0; r < UHD; r ++) {
                NHDS[r] = 0.0;
                NCSP[r] = 0.0;
                nodeOfUHD[r] = 0;
        }

        for(int n = 0; n < noinw; n ++) {
                if(relRC[n] != 0) {
                        nodeOfUHD[k] = n + 1;
                        NHDS[k] = HDS[n];
                        NCSP[k] = CSP[n];
                        k = k + 1;
                }
        }       
}

void solveSimulEqns(double coefs[], double x[], double c[], const int UHD) {
        //Gauss-Jordan Reduction Technique              
        int ind, id, nd, id1, id2, id3;
        for(int p = 0; p < UHD; p++)
                x[p] = 0.0;
        for(int i = 0; i < UHD; i++) {
                nd = UHD*i + i;         
                if(coefs[nd] != 0.0) {
                        for(int h = 0;h < UHD; h++) {
                                ind = UHD*i + h;        
                                coefs[ind] = coefs[ind]/coefs[nd];
                        }
                        c[i] = c[i]/coefs[nd];
                }

                for(int j = 0; j < UHD; j++) {
                        id1 = UHD*j + i;                        
                        if(coefs[id1] < 0.0)
                                coefs[id1] = -(coefs[id1]);
                        if(j != i) {
                                id2 = UHD*j + i;        
                                if(coefs[id2] < 0.0) {
                                       for(int l = i; l < UHD; l++) {
                                                id3 = UHD*j + l;
                                                id = UHD*i + l;         
                                                coefs[id3] = coefs[id3] + coefs[id]*coefs[id1];
                                        }
                                        c[j] = c[j] + c[i]*coefs[id1];
                                } else if(coefs[id2] > 0.0) {
                                        for(int l = i; l < UHD; l++) {
                                                id3 = UHD*j + l;
                                                id = UHD*i + l; 
                                                coefs[id3] = coefs[id3] - coefs[id]*coefs[id1];
                                        }
                                        c[j] = c[j] - c[i]*coefs[id1];
                                }
                        }
                }
        }
}

void initpdRead(double pdRead[],const int nenw)  {
        for(int k = 0; k < nenw; k++)
                pdRead[k] = 0.0;
}

void insertTOpdRead(double pdRead[], const int elem, const double diam)  {
        if(pdRead[(elem - 1)] == 0.0)
                pdRead[(elem - 1)] = diam;
        
}

void initplRead(double plRead[],const int nenw)  {
        for(int k = 0; k < nenw; k++)
                plRead[k] = 0.0;
}


void insertTOplRead(double plRead[], const int elem, const double plen)  {
        if(plRead[(elem - 1)] == 0.0)
                plRead[(elem - 1)] = plen;
        
}


void readInt(int& inv) {
        char s[15];        
        int thro; 
        do {
                
                std::cin >> s;
                inv = atoi(s);             

                std::cout << std::endl << "     Data entered is " << inv;
                std::cout << std::endl;
                std::cout << std::endl << "     Is entry correct (Y/N)? ";

                char ans;               
                std::cin >> ans;
                if((ans == 'N') || (ans == 'n'))                    
                	thro = 1;                       
                else if ((ans == 'Y') || (ans == 'y')) 
                        thro = 0;
                else
                        thro = 1;               

                if (thro > 0) 
                        std::cout << std::endl << "     ERROR: Re-enter data: ";
        } while (thro > 0);
}


void readNode(int& inv, const int noinw) {
        char s[15];        
        int thro, flag; 
        do {
                
                std::cin >> s;
                inv = atoi(s);
                if ((inv <= 0) || (inv > noinw)) 
                    flag = 1;
                else if ((inv > 0) && (inv <= noinw))
                    flag = 0;
                else 
                    flag = 1;                    

                std::cout << std::endl << "     Data entered is " << inv;
                std::cout << std::endl;
                std::cout << std::endl << "     Is entry correct (Y/N)? ";

                char ans;               
                std::cin >> ans;
                if((ans == 'N') || (ans == 'n'))                    
                        thro = 1;                       
                else if ((ans == 'Y') || (ans == 'y')) 
                        thro = 0;
                else
                        thro = 1;               

                if ((flag > 0) || (thro > 0)) 
                        std::cout << std::endl << "     ERROR: Re-enter data: ";
        } while ((flag > 0) || (thro > 0));
}


void readElem(int& inv, const int nenw) {
        char s[15];        
        int thro, flag; 
        do {
                
                std::cin >> s;
                inv = atoi(s);
                if ((inv <= 0) || (inv > nenw)) 
                    flag = 1;
                else if ((inv > 0) && (inv <= nenw))
                    flag = 0;
                else 
                    flag = 1;                    

                std::cout << std::endl << "     Data entered is " << inv;
                std::cout << std::endl;
                std::cout << std::endl << "     Is entry correct (Y/N)? ";

                char ans;               
                std::cin >> ans;
                if((ans == 'N') || (ans == 'n'))                    
                        thro = 1;                       
                else if ((ans == 'Y') || (ans == 'y')) 
                        thro = 0;
                else
                        thro = 1;               

                if ((flag > 0) || (thro > 0)) 
                        std::cout << std::endl << "     ERROR: Re-enter data: ";
        } while ((flag > 0) || (thro > 0));
}


void readDouble(double& flov) {
        char s[15];        
        int thro; 	//flag;
        do {
                
                std::cin >> s;
                flov = atof(s);
                std::cout << std::endl << "     Data entered is " << flov;
                std::cout << std::endl;
                std::cout << std::endl << "     Is entry correct (Y/N)? ";

                char ans;               
                std::cin >> ans;
                if((ans == 'N') || (ans == 'n'))                    
                        thro = 1;                       
                else if ((ans == 'Y') || (ans == 'y')) 
                        thro = 0;
                else
                        thro = 1;                             

		if (thro > 0) 
                        std::cout << std::endl << "     ERROR: Re-enter data: ";
        } while (thro > 0);
}

void readCHD(double& flov) {
        char s[15];        
        int thro; 
        do {                
                std::cin >> s;
                flov = atof(s);                
                std::cout << std::endl << "     Data entered is " << flov;
                std::cout << std::endl;
                std::cout << std::endl << "     Is entry correct (Y/N)? ";

                char ans;               
                std::cin >> ans;
                if((ans == 'N') || (ans == 'n'))                    
                        thro = 1;                       
                else if ((ans == 'Y') || (ans == 'y')) 
                        thro = 0;
                else
                        thro = 1;               

                if (thro > 0) 
                        std::cout << std::endl << "     ERROR: Re-enter data: ";
        } while (thro > 0);
}

void NHIncrement(double H[], double HP[],const double NCSP[], const int nodeOfUHD[],const int UHD, const int noinw) {
        int j = 0;      
        int v = 0;

        for(int c = 0; c < noinw; c++)          
                        HP[c] = H[c];
                
        for(int n = 0; n < noinw; n++) {
                if (v < UHD) {
                        j = nodeOfUHD[v];
                        j = j - 1;
                        H[j] = H[j] + NCSP[v];
                } else if (v == UHD) 
                        break;
                v = v + 1;
        }                              
}

void ResInConsum(double CSP2[], const double CSP3[],const int noinw) {
        for(int n = 0; n < noinw; n++)
                CSP2[n] = CSP3[n];
        
}

void ResInHead(double HDS[], const double HD[], const int noinw) {
        for(int n = 0; n < noinw; n++)                  
                HDS[n] = HD[n];
        
}

void multi(const double HDS[], double CSP[], const double KVals[],const int noinw) {
        int ind;
        for(int v = 0; v < noinw; v++)
                CSP[v] = 0.00;
        for(int i = 0; i < noinw; i++) {
                for(int j = 0; j < noinw; j++) {
                        ind = noinw*i + j;                                                              
                        CSP[i] = CSP[i] + KVals[ind]*HDS[j];
                }
        }
}

void Results(const double kel[], double NDS[], double MDS[], const double H[],const int nenw,const int noinw)  {
        int nodnu = 0;
        int nuelem = 0;
        int bnode = 0;
        int enode = 0;
        int elem = 0;
        double plen = 0.0;
        double pdiam = 0.0;
        double ffact = 0.0;
        double a, hb, he;
        int b, e;
        std::ifstream inf2("zon12.dat", std::ios::in);
        std::ifstream inf3("zon13.dat", std::ios::in);
        for(int k = 0; k < noinw; k++) 
                NDS[k] = 0.0;
        for(int l = 0; l < nenw; l++) 
                MDS[l] = 0.0;
                
        int j = 0;      
        while(inf2 >> nodnu >> nuelem) {
                j = nodnu;
                int cnt = 0;
                double sumMDS = 0.0;
                int nk = 0.0;
                while(inf3 >> elem >> bnode >> enode >>  plen >> pdiam >> ffact) {
                        nk = elem;                      
                        cnt = cnt + 1;                  
                        b = bnode; // Beginning Node No.
                        e = enode; // Ending Node No.
                        hb = H[(b - 1)];  // Head of beginning Node
                        he = H[(e - 1)];  // Head of ending Node                        
                        a = hb - he;                    
                        double sgn;
                        if (a >= 0.0)
                                sgn = 1.0;
                        else if (a < 0.0)                               
                                sgn = -1.0;
                        
                        MDS[nk - 1] = sgn*kel[nk - 1]*(pow(fabs(a), 0.50));
                        sumMDS = sumMDS + MDS[nk - 1];                  
                        if (cnt == nuelem)
                                break;

                }
                NDS[j - 1] = sumMDS;
                j = j + 1;
        }       
        inf2.close();
        inf3.close();
}

void OutPut(const double H[], const double MDS[], const double NDS[],const int nenw, const int noinw)  {
        std::cout << "   " << std::endl;  
        std::cout << std::endl;
        int n = 0; 
        std::cout << "     Node Number     " <<  "  Nodal Heads       " << std::endl;

        std::ofstream outf11("zon111.dat", std::ios::out);        
        for(n = 0; n < noinw; n++) {    
                std::cout << "         " << n + 1 <<  "                " << H[n] << std::endl;
                outf11 << n + 1 <<  " " << H[n] << std::endl;
                //std::cin >> xx;                            
        }
        outf11.close();
        std::cout << std::endl;        
    
        int c;
	std::cout << "     Press any key to continue  " << std::endl;
	c = getchar();

        std::ofstream outf12("zon112.dat", std::ios::out);
        std::cout << "     Node Number     " <<  "  Nodal Discharges  " <<  std::endl;
        for(n = 0; n < noinw; n++) {    
                std::cout << "         " << n + 1 <<  "                " << NDS[n] << std::endl;
                outf12 << n + 1 <<  " " << NDS[n] << std::endl;
                //std::cin >> xx;            
        }       
        outf12.close();
        std::cout << std::endl;

        std::cout << "     Press any key to continue  " << std::endl;
        c = getchar();
        

        std::ofstream outf13("zon113.dat", std::ios::out);
        std::cout << "     Pipe Number  " <<  "  Pipe Discharges " << std::endl;
        for(n = 0; n < nenw; n++) {     
                std::cout << "         " << n + 1 <<  "                " << MDS[n] << std::endl;
                outf13 << n + 1 <<  " " << MDS[n] << std::endl;
                //std::cin >> xx;            
        }               
        outf13.close();
        std::cout << std::endl;
        
	std::cout << "     Press any key to continue  " << std::endl;
        c = getchar();     
}

void AllNodalHeads(double H[], const double HD[],const int noinw) {
        for(int n = 0; n < noinw; n++){ 
                if (HD[n] > 0.0)
                        H[n] = HD[n];
        }
}

void CopyOfKVals(double DKVals[], const double KVals[], const int UHD, const int noinw) {
        int ind;
        for(int i = 0; i < noinw; i++) {                                
                for(int j = 0; j < noinw; j++) {
                        ind = UHD*i + j;
                        DKVals[ind] = KVals[ind];
                }
        }
}

void ArraysSet1 (double* &ncoef,double* &CSP, double* &HDS,double* &HD,
         double* &CSP2, double* &CSP3,double* &NDS,double* &H,double* &HP,const int noinw) {
        ncoef = new double[noinw];
        CSP = new double[noinw];
        HDS = new double[noinw];
        HD = new double[noinw];
        CSP2 = new double[noinw];
        CSP3 = new double[noinw];
        NDS = new double[noinw];
        H = new double[noinw];
        HP = new double[noinw];
        //x = new double[noinw];
        //y = new double[noinw];  
}

void ArraysSet2 (double* &KVals,double* &DKVals, int* &ndnu, int* &relRC, const int noinw) {
        int nm;
        nm = noinw*noinw;       
        KVals = new double[nm];
        DKVals = new double[nm];
        ndnu = new int[noinw];
        relRC = new int[noinw]; 
}


void ArraysSet3 (double* &plRead,double* &pdRead,double* &MDS, double* &kel, double* &epl, double* &epd, double* &eff,
         int* &elnu, int* &bno,int* &eno, const int nenw) {
        plRead = new double[nenw];
        pdRead = new double[nenw];
        MDS = new double[nenw];
        kel = new double[nenw];
        epl = new double[nenw];
        epd = new double[nenw];
        eff = new double[nenw];
        elnu = new int[nenw];
        bno = new int[nenw];
        eno = new int[nenw];            
}

void ArraysSet4 (double* &NHDS,double* &NCSP, double* &NKVals, int* &nodeOfKHD, 
        int* &nodeOfUHD, const int UHD, const int NPH) {
        int vno;
        NHDS = new double[UHD];
        NCSP = new double[UHD];
        vno = UHD*UHD; 
        NKVals = new double[vno];
        nodeOfKHD = new int[NPH];
        nodeOfUHD = new int[UHD];
}

void NetWorkAttributes (int& noinw, int& nenw) {
        int nnode, nelem;   
        std::ifstream inf5("zon15.dat", std::ios::in);
        if(!inf5) {                
                ReadNWAttributes (noinw, nenw);
        } else {
                inf5 >> nnode >> nelem;         
                inf5.close();           
                noinw = nnode;
                nenw = nelem;
        }
        
}


void ReadIteTol(double& tol, int& mnite) {
        std::cout << std::endl;
        std::cout << std::endl;
        std::cout << "     ITERATION AND CONVERGENCE PARAMETERS" << std::endl;
        std::cout << std::endl;                   
        std::cout << "     Enter number of iterations:  ";
        //std::cin >> mnite;
        readInt(mnite);
        std::cout << std::endl;
        std::cout << "     Enter tolerance value:  ";
        readDouble(tol);
}


void ReadNWAttributes (int& noinw, int& nenw) {
        int nnode, nelem;     //, NKHD, NCOS, NUHD;     
        std::ofstream outf5("zon15.dat", std::ios::out);  
        
        std::cout << "     Enter number of nodes in the network: ";
        //std::cin >> nnode;
        readInt(nnode);
        noinw = nnode;
        std::cout << std::endl;
        std::cout << "     Enter number of elements in the network: ";
        readInt(nelem);
        nenw = nelem;
        std::cout << std::endl;
        outf5 << nnode <<  " " << nelem << std::endl;  
        outf5.close();
}


void outputData(const int ndnu[], const double HDS[], const double CSP[], const int elnu[], const int bno[],
        const int eno[], const double epl[], const double epd[], const double eff[], const int nenw, const int noinw) {
        int c;
        std::cout << std::endl;
        int n = 0;

        n = 0;
        std::cout << std::endl;   
        std::cout << "     Pipe Number  " << "Start Node  " <<  "End Node  " << "P. Length  "  << "P. Diameter  " <<  " F/Factor  " << std::endl;
        std::ofstream outf9("zon19.dat", std::ios::out);
        for(n = 0; n < nenw; n++) { 
                std::cout << "        " << elnu[n]  << "                "<< eno[n]  <<  "          "<< bno[n]  << "       " << epl[n] <<  "        " << epd[n] <<  "        "<< eff[n] << std::endl;
                outf9 << elnu[n] <<  " " << eno[n] << " " << bno[n] << " " << epl[n] << " " << epd[n] << " " << eff[n]  << std::endl;
        }
        outf9.close();
        std::cout << std::endl;
        
	std::cout << "     Press any key to continue  " << std::endl;
        c = getchar();     
}


void PrintHeader() {
    char str0[] = "               ";
    fprintf(stdout, "%s\n", str0);
    fprintf(stdout, "%s\n", str0);
    char str[] = "               NETWORK ANALYSIS BY NON-LINEAR METHOD USING";
    fprintf(stdout, "%s\n", str);
    char str1[] = "                              DARCY-WEISBACH EQUATION";
    fprintf(stdout, "%s\n", str1);
    char str2[] = "          ";
    fprintf(stdout, "%s\n", str2);
    fprintf(stdout, "%s\n", str2);
}

void PrintIteTol(const double tol, const int mnite) {
    char str[] = "     ITERATION AND CONVERGENCE PARAMETERS";
    fprintf(stdout, "%s\n", str);
    char str1[] = "     Number of iterations:              ";
    fprintf(stdout, "%s %9d\n", str1, mnite);
    char str2[] = "     Tolerance:                         ";
    fprintf(stdout, "%s %9.4f\n", str2, tol);
}

void PrintErrConv(const double error, const int istep, const double tol) {
    char str[] = "     Error:                              ";
    fprintf(stdout, "%s %8.4f\n", str, error);    
    if ((error - tol) <= 0.0) {
        char str1[] = "     Convergence achieved on cycle: ";
        fprintf(stdout, "%s %13d\n", str1, istep);       
    } else if ((error - tol) > 0.0) {
        char str2[] = "     Maximum number of iteration reached without convergence";        
        fprintf(stdout, "%s\n", str2);
    } 
    char str3[] = "      ";
    fprintf(stdout, "%s\n", str3);
}

void PrintData(const int ndnu[], const double HDS[], const double CSP[], const int elnu[], const int bno[],
        const int eno[], const double epl[], const double epd[], const double eff[], const int nenw, const int noinw) {
        char str1[] = "      ";
        fprintf(stdout, "%s\n", str1);
        int n = 0;
        int nn = 0;
	/*
        char str[] = "     Node Number     Head         Consumption  ";   
	fprintf(stdout, "%s\n", str);   
        for(int v = 0; v < noinw; v++) 
            fprintf(stdout, "%10d %16.4f %16.4f\n", ndnu[v], HDS[v], CSP[v]);
            
        fprintf(stdout, "%s\n", str1);
	*/
        char str2[] = "     Pipe Number  Start Node  End Node  P. Length  P. Diameter   F/Factor  ";  
	fprintf(stdout, "%s\n", str2);      
        for(n = 0; n < nenw; n++) 
            fprintf(stdout, "%10d %15d %10d %14.4f %12.4f %12.4f\n", elnu[n], eno[n], bno[n], epl[n], epd[n], eff[n]);
        
        fprintf(stdout, "%s\n", str1);    
}


void PrintResult(const double H[], const double MDS[], const double NDS[],const int nenw, const int noinw)  {
        char str1[] = "      ";
        fprintf(stdout, "%s\n", str1);
        int n = 0;      
        int nn = 0;
        char str2[] = "     Node Number       Nodal Heads       ";
	fprintf(stdout, "%s\n", str2);              
        for(n = 0; n < noinw; n++) {
            nn = n + 1;  
            fprintf(stdout, "%10d %18.4f\n", nn, H[n]);                                           
        }        
        fprintf(stdout, "%s\n", str1);
        char str3[] = "     Node Number       Nodal Discharges  ";
	fprintf(stdout, "%s\n", str3);
        for(n = 0; n < noinw; n++) {
            nn = n + 1;
            fprintf(stdout, "%10d %18.4f\n", nn, NDS[n]);                          
        }       
        fprintf(stdout, "%s\n", str1);
        char str4[] = "     Pipe Number    Pipe Discharges ";
	fprintf(stdout, "%s\n", str4);
        for(n = 0; n < nenw; n++) {
            nn = n + 1; 
            fprintf(stdout, "%10d %18.4f\n", nn, MDS[n]);           
        }               
        fprintf(stdout, "%s\n", str1);    
}

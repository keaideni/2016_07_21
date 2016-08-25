#include "QWave.h"
#include "Corr.h"
#include "DMRGP.h"
//================function define===================================================
void CalcuCorr(const int& OrbitalM, const QWave& fwave, std::ofstream& Fdata);
double CacuCorr(const OP& corrn, const OP& corrc, const OP& corrcdag, const QWave& fwave);
double CacuCorr(const OP& corrc, const OP& corrcdag, const QWave& fwave);
double CacuCorrM(const OP& corrc, const OP& corrcdag, const QWave& fwave);


//=================calculate the product of two wave====================================
double prod2wave(const QWave& wave1, const QWave& wave2);
double prod2wave(const QWave& wave1, const QWave& wave2)
{
        double prod;
        for(auto it = wave1.WavePart.begin(); it != wave1.WavePart.end(); ++it)
        {
                auto itt = wave2.WavePart.find(it->first);
                if(itt != wave2.WavePart.end())
                {
                        for(auto matit1 = it->second.QMat.begin(); matit1 != it->second.QMat.end(); ++matit1)
                        {
                                auto matit2 = itt->second.QMat.find(matit1->first);
                                if(matit2 != itt->second.QMat.end())
                                {
                                        if(it->second.RLQ.at(matit1->first) == itt->second.RLQ.at(matit1->first))
                                        {
                                                for(int i = 0; i < matit1->second.rows(); ++i)
                                                {
                                                        for(int j = 0; j < matit1->second.cols(); ++j)
                                                        {
                                                                prod += matit1->second(i,j)*matit2->second(i,j);
                                                        }
                                                }
                                        }
                                }
                        }
                }
        }
        return prod;
}
//================================================================================================================


//=================calculate two types of correlation function========================
void CalcuCorr(const int& OrbitalM, const QWave& fwave, std::ofstream& Fdata)
{
        std::ofstream outfile1, outfile11;//to save the normalized correlation.
        std::ofstream outfile2, outfile22;//to save the not normalized correlation.
        if(OrbitalM % 2 == 1)
        {
                                        
                outfile1.open("./result/resonatorN.txt");//to store the type of particle besides the OrbitalM.
                outfile2.open("./result/resonator.txt");//to store the type of particle besides the OrbitalM.
                
                outfile11.open("./result/qubitN.txt");//to store the type of particle next besides the OrbitalM.
                outfile22.open("./result/qubit.txt");//to store the type of particles next besides the orbitalM.
        }else
        {
                                        
                outfile1.open("./result/qubitN.txt");//to store the type of particle besides the OrbitalM.
                outfile2.open("./result/qubit.txt");//to store the type of particles besides the orbitalM.
                
                outfile11.open("./result/resonatorN.txt");//to store the type of particle next besides the OrbitalM.
                outfile22.open("./result/resonator.txt");//to store the type of particles next besides the orbitalM.
        }



        int i(OrbitalM + 1);//this label for Sigma;
        int j(OrbitalM - 1);//this label for Sigmadag and N;
        int fflag(1);

        double CorrLenth(0);
        double CorrSum(0);
        double correlation(0);

        Corr corr, corrdag, corrn;
        while(true)
        {
                                        
                                        

                corr.read(i, 1);corrdag.read(j, 2);corrn.read(j, 3);
                                        //corrn.show();Sys.SubSysEye.show();
                //corr.show();
                correlation = CacuCorr(corrn.CorrO, corr.CorrO, corrdag.CorrO, fwave);


                int distence(i - j);

                /*std::cout<<"Distence = " << distence << ", the correlation = "
                <<correlation<<std::endl;*/
                outfile1<<"Distence = " << distence << ", the correlation = "
                <<correlation<<std::endl;


                correlation = CacuCorr(corr.CorrO, corrdag.CorrO, fwave);
                /*std::cout<<"Distence = " << distence << ", the correlation = "
                <<correlation<<std::endl;*/
                outfile2<<"Distence = " << distence << ", the correlation = "
                <<correlation<<std::endl;

                CorrLenth += pow((i-j)/2, 2)*correlation;
                
                CorrSum +=correlation;
        

                if(fflag == 1)
                {
                        i += 2;
                }else
                {
                        j-= 2;
                }
                if((i>OrbitalM* 2)||(j<=0))break;
                fflag *= -1;


        }




        
        CorrLenth = CorrLenth/CorrSum;
        if(OrbitalM %2 == 0)
        {
                Fdata << "the Qubit correlation is " <<CorrLenth<<std::endl;
                //std::cout << "the Qubit correlated length = " <<CorrLenth<<std::endl;
        }else
        {
                Fdata << "the Resonator correlation is " <<CorrLenth<<std::endl;
                //std::cout << "the Resonator correlated length = " <<CorrLenth<<std::endl;
        }



        i=(OrbitalM + 2);//this label for Sigma;
        j=(OrbitalM - 2);//this label for Sigmadag and N;
        fflag=1;
        CorrLenth = 0;
        CorrSum = 0;
        while(true)
        {
                                        
                                        

                corr.read(i, 1);corrdag.read(j, 2);corrn.read(j, 3);
                //corrn.show();Sys.SubSysEye.show();

                correlation = CacuCorr(corrn.CorrO, corr.CorrO, corrdag.CorrO, fwave);

                int distence(i - j);

                /*std::cout<<"Distence = " << distence << ", the correlation = "
                <<correlation<<std::endl;*/
                outfile11<<"Distence = " << distence << ", the correlation = "
                <<correlation<<std::endl;


                correlation = CacuCorr(corr.CorrO, corrdag.CorrO, fwave);
                /*std::cout<<"Distence = " << distence << ", the correlation = "
                <<correlation<<std::endl;*/
                outfile22<<"Distence = " << distence << ", the correlation = "
                <<correlation<<std::endl;
                CorrLenth += pow((i-j)/2, 2)*correlation;
                CorrSum += correlation;

                if(fflag == 1)
                {
                        i += 2;
                }else
                {
                        j-= 2;
                }
                if((i>OrbitalM*2)||(j<=0))break;
                fflag *= -1;


        }

        //============to add the site OrbitalM============================
        Parameter para;
        para.read();
        Sub m(para, OrbitalM);
        corr.read(OrbitalM+2, 1);

        correlation=CacuCorrM(corr.CorrO, m.SubSysCdag, fwave);

        
        CorrLenth += correlation;
                
        CorrSum +=correlation;
        //================================================================

        CorrLenth /= CorrSum;

        if(OrbitalM %2 == 0)
        {
                Fdata << "the Resonator correlation = " <<CorrLenth<<std::endl;
                //std::cout << "the Resonator correlated length = " <<CorrLenth<<std::endl;
        }else
        {
                Fdata << "the Qubit correlation = " <<CorrLenth<<std::endl;
                //std::cout << "the Qubit correlated length = " <<CorrLenth<<std::endl;
        }

        outfile1.close();
        outfile2.close();
        outfile11.close();
        outfile22.close();
}
//=================================================================================================================


//=======================the first correlation which is normalized===========================================
double CacuCorr(const OP& corrn, const OP& corrc, const OP& corrcdag, const QWave& fwave)
{
        

        QWave wave1, wave2;
        //std::vector<double> f, f1;
        //fwave.Wave2f(f);
        //QWave ffwave(fwave);
        //double number(0);
        //double correlation(0);
        wave1.clear();
        wave1.OSWave2New(corrn, fwave);
        //ffwave.initial(wave1);
        //ffwave.Wave2f(f1);

        double number(prod2wave(fwave, wave1));

        /*for(int i = 0; i < f.size(); ++i)
        {
                number += f[i]*f1[i];
        }*/

        wave1.clear();
        wave1.OSWave2New(corrcdag, fwave);
        wave2.clear();
        wave2.OEWave2New(corrc, wave1);
        //ffwave.initial(wave2);
        //ffwave.Wave2f(f1);
        /*correlation = 0;
        for(int i = 0; i < f.size(); ++i)
        {
                correlation += f[i]*f1[i];
        }
        correlation /= number;*/
        //std::cout<<correlation<<std::endl;

        double correlation(prod2wave(fwave, wave2)/number);
        
        return correlation;
        //std::cout<<number<<std::endl;


}
//=====================================================================================================

//================the second correlation which isn't normalized========================================
double CacuCorr(const OP& corrc, const OP& corrcdag, const QWave& fwave)
{
        QWave wave1, wave2;
        //std::vector<double> f, f1;
        //fwave.Wave2f(f);
        //QWave ffwave(fwave);
        
        //double correlation(0);

        wave1.clear();
        wave1.OSWave2New(corrcdag, fwave);
        wave2.clear();
        wave2.OEWave2New(corrc, wave1);
        /*ffwave.initial(wave2);
        ffwave.Wave2f(f1);
        correlation = 0;
        for(int i = 0; i < f.size(); ++i)
        {
                correlation += f[i]*f1[i];
        }*/

        return prod2wave(fwave, wave2);
        
        
        //return correlation;
        //std::cout<<number<<std::endl;


}
//====================used to calculate the correlation length which is defined directly================



//===============to add the OrbitalM site for the correlation at length 2================================
double CacuCorrM(const OP& corrc, const OP& corrcdag, const QWave& fwave)
{
        QWave wave1, wave2;
        //std::vector<double> f, f1;
        //fwave.Wave2f(f);
        //QWave ffwave(fwave);
        
        //double correlation(0);

        wave1.clear();
        wave1.OMWave2New(corrcdag, fwave);
        wave2.clear();
        wave2.OEWave2New(corrc, wave1);
        /*ffwave.initial(wave2);
        ffwave.Wave2f(f1);
        correlation = 0;
        for(int i = 0; i < f.size(); ++i)
        {
                correlation += f[i]*f1[i];
        }*/

        return prod2wave(fwave, wave2);
        
        
        //return correlation;
        //std::cout<<number<<std::endl;


}
//============================================================================================



//===============calculate the density========================================================
double density(const OP& corrn, const QWave& fwave, const int& OrbitalM, const int& orbital);
double density(const OP& corrn, const QWave& fwave, const int& OrbitalM, const int& orbital)
{

        //double density(0);
        //QWave ffwave(fwave);
        QWave wave1;
        //std::vector<double> f1, f2;

        //fwave.Wave2f(f1);

        if(orbital<OrbitalM)
        {
                wave1.OSWave2New(corrn, fwave);
        }else if(orbital > OrbitalM)
        {
                wave1.OEWave2New(corrn, fwave);
        }else
        {
                wave1.OMWave2New(corrn, fwave);
        }

        /*ffwave.initial(wave1);
        ffwave.Wave2f(f2);

        for(int i = 0; i < f2.size(); ++i)
        {
                density += f2.at(i)*f1.at(i);
        }*/

        return prod2wave(fwave, wave1);

        //return density;

}

void calcudensity(const int& OrbitalM, const QWave& fwave, const int& endN);
void calcudensity(const int& OrbitalM, const QWave& fwave, const int& endN)
{

        std::ofstream Qdensity("./result/Qdensity");
        std::ofstream Rdensity("./result/Rdensity");
        Corr corrn;
        for(int i = 1; i < OrbitalM; ++i)
        {
                corrn.read(i, 3);
                double den(density(corrn.CorrO, fwave, OrbitalM, corrn.orbital()));
                if(i%2 == 1)
                {
                        Qdensity<<"site = " << i << ", the density = "
                        <<den<<std::endl;
                        /*std::cout<<"site = " << i << ", the density = "
                        <<den<<std::endl;*/
                }else
                {
                        Rdensity<<"site = " << i << ", the density = "
                        <<den<<std::endl;
                        /*std::cout<<"site = " << i << ", the density = "
                        <<den<<std::endl;*/
                }

        }
        Parameter para;
        para.read();
        Sub m(para, OrbitalM);

        OP temp;
        temp.time(m.SubSysCdag, m.SubSysC);

        double den(density(temp, fwave, OrbitalM, OrbitalM));
        if(OrbitalM%2 == 1)
        {
                Qdensity<<"site = " << OrbitalM << ", the density = "
                <<den<<std::endl;
                /*std::cout<<"site = " << OrbitalM << ", the density = "
                <<den<<std::endl;*/
        }else
        {
                Rdensity<<"site = " << OrbitalM << ", the density = "
                <<den<<std::endl;
                /*std::cout<<"site = " << OrbitalM << ", the density = "
                <<den<<std::endl;*/
        }


        for(int i = OrbitalM+1; i <= endN; ++i)
        {
                corrn.read(i, 3);
                double den(density(corrn.CorrO, fwave, OrbitalM, corrn.orbital()));
                if(i%2 == 1)
                {
                        Qdensity<<"site = " << i << ", the density = "
                        <<den<<std::endl;
                        /*std::cout<<"site = " << i << ", the density = "
                        <<den<<std::endl;*/
                }else
                {
                        Rdensity<<"site = " << i << ", the density = "
                        <<den<<std::endl;
                        /*std::cout<<"site = " << i << ", the density = "
                        <<den<<std::endl;*/
                }

        }

        Qdensity.close();
        Rdensity.close();


}

//===================================================================================================



//====================calculate the structure factor==================================================
double QdenCorr(const QWave& fwave, const int& OrbitalM, const int& orbital);
double QdenCorr(const QWave& fwave, const int& OrbitalM, const int& orbital)
{
        QWave wave1;
        if(orbital < OrbitalM)
        {
                Corr corr;
                corr.read(orbital, 5);
                wave1.OSWave2New(corr.CorrO, fwave);
        }else if(orbital > OrbitalM)
        {
                Corr corr1, corr2;
                corr1.read(1, 3);corr2.read(orbital, 3);
                QWave wave2;
                wave2.OEWave2New(corr2.CorrO, fwave);
                wave1.OSWave2New(corr1.CorrO, wave2);
        }else
        {
                Parameter para;
                para.read();

                Sub m(para, OrbitalM);

                OP tempm;
                tempm.time(m.SubSysCdag, m.SubSysC);

                QWave wave2;
                wave2.OMWave2New(tempm, fwave);

                Corr corr1;
                corr1.read(1, 3);
                wave1.OSWave2New(corr1.CorrO, wave2);
        }


        return prod2wave(wave1, fwave);
}



double RdenCorr(const QWave& fwave, const int& OrbitalM, const int& orbital);
double RdenCorr(const QWave& fwave, const int& OrbitalM, const int& orbital)
{
        QWave wave1;
        if(orbital < OrbitalM)
        {
                Corr corr;
                corr.read(orbital, 7);
                wave1.OSWave2New(corr.CorrO, fwave);
        }else if(orbital > OrbitalM)
        {
                Corr corr1, corr2;
                corr1.read(2, 3);corr2.read(orbital, 3);
                QWave wave2;
                wave2.OEWave2New(corr2.CorrO, fwave);
                wave1.OSWave2New(corr1.CorrO, wave2);
        }else
        {
                Parameter para;
                para.read();

                Sub m(para, OrbitalM);

                OP tempm;
                tempm.time(m.SubSysCdag, m.SubSysC);

                QWave wave2;
                wave2.OMWave2New(tempm, fwave);

                Corr corr1;
                corr1.read(2, 3);
                wave1.OSWave2New(corr1.CorrO, wave2);
        }


        return prod2wave(wave1, fwave);
}



double calcustructure(const QWave& fwave, const int& OrbitalM, const int& endN, std::ofstream& Fdata);
double calcustructure(const QWave& fwave, const int& OrbitalM, const int& endN, std::ofstream& Fdata)
{
        double Sq(QdenCorr(fwave, OrbitalM, 1)), Sr(RdenCorr(fwave, OrbitalM, 2));
        int flag(-1), fflag(-1);



        for(int i=3; i<=endN; ++i)
        {
                if(i%2 == 1)
                {
                        Sq += flag*2*QdenCorr(fwave, OrbitalM, i);
                        flag*=-1;
                }else
                {
                        Sr += fflag*2*RdenCorr(fwave, OrbitalM, i);
                        fflag*=-1;
                }
        }
        Sq += flag*QdenCorr(fwave, OrbitalM, endN+1);
        Sr += fflag*RdenCorr(fwave, OrbitalM, endN+2);


        Sq /= endN;
        Sr /= endN;

        Fdata<<"the Qubit structure factor is "<<Sq<<std::endl;
        Fdata<<"the Resonator structure factor is "<<Sr<<std::endl;
}
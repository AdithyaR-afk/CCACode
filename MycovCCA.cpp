// MycovCCA.cpp : Defines the entry point for the console application.
//




#include "stdafx.h"
#include<conio.h>
#include<iostream>
#include<cmath>
#include <Eigen/Dense>


#include <Eigen/Eigenvalues>


using namespace std;
using namespace Eigen;
using Eigen::MatrixXd;
int main()
{
        int szar,szac,szbr,szbc;
        //Enter the size of the rows and columns of matrices a and b.
        szar=5;
        szac=5;
        szbr=5;
        szbc=5;
        MatrixXf a(szar,szac);
        a<< 5,9,4,1,2,
                7,9,7,2,4,
                7,9,7,3,6,
                5,6,7,8,8,
                6,5,4,3,4;
        MatrixXf b(szbr,szbc);
        b<< 6,9,4,4,8,
                4,4,8,4,7,
                4,4,8,4,6,
                8,9,2,3,5,
                5,4,5,6,7;


        
        //Matrix c is the concatenated version of a and b appended from the bottom
        MatrixXf c(szar+szbr,szac);
        //This for loop block does the concatenation.
     for(int s=0; s<(szar+szbr);s++)
        {
                if(s<szar){
                for(int j=0; j<szac;j++)
                {
                        c(s,j)=a(s,j);
                        //cout<<c(s,j)<< "a"<<endl;
                }}
                        if(s>szar-1){
                for(int j=0; j<szac;j++)
                {
                        c(s,j)=b(s-szbr,j);
                        //cout<<c(s,j)<< "b"<<endl;
                        


                }}}cout<<endl<<endl<<"concatenated matrix"<<endl<<c<<endl<<endl;
        int rowsum= szar+szbr;
        //Initialize a zero array v
        //ArrayXf k = ArrayXf::Zero(rowsum);
        ArrayXf v = ArrayXf::Zero(rowsum);
        int nrowsc= c.rows();
        int nrowsa=a.rows();
        int ncolsa=a.cols();
        //This block computes the rowwise mean of each row in matrix c.
        for(int i=0; i<(rowsum);i++)
        {
                for(int w=0;w<szac;w++)
                {
        
          v[i]=v[i]+c(i,w);


                }
                
                v[i]=v[i]/ncolsa;
                
                
        }
        cout<<"mean array "<<endl<<v<<endl<<endl;
        //Initializing the covariance matrix 'cca' with zeroes.
        MatrixXf cca(2*nrowsa,2*nrowsa);
        for(int p=0;p<(2*nrowsa);p++)
        {
                for(int q=0; q<(2*nrowsa);q++)
                {
                        cca(p,q)=0;
                }
        }


        cout<<"rowsum: "<<rowsum<<endl<<endl<<"cca"<<endl<<cca<<endl;
        //Computing the values of the covariance matrix.
        for(int x=0; x<rowsum;x++)
        {
                for(int y=0; y<(2*nrowsa);y++)
                {
                        for(int z=0; z<ncolsa;z++)
                        {
                                cca(x,y) = cca(x,y) + (c(x,z)-v[x])*(c(y,z)-v[y]);
                        
                        }
                        cca(x,y)=cca(x,y)/(ncolsa-1);
                        
                }
        }


        cout<<endl<<" covariance matrix 'cca' output"<<endl<<cca<<endl;
        
//Extracting each matrix from covariance matrix
        MatrixXf cxx = cca.block(0,0,nrowsa,nrowsa);
        MatrixXf cxy = cca.block(0,nrowsa,nrowsa,nrowsa);
        MatrixXf cyx = cca.block(nrowsa,0,nrowsa,nrowsa);
        MatrixXf cyy = cca.block(nrowsa,nrowsa,nrowsa,nrowsa);


        //Generating an identity matrix.
        MatrixXf id = MatrixXf::Identity(nrowsa,nrowsa);
        
        //In order to compute the inverse for cxx and cyy matrices, a small value is added to make sure that the inverse exists.
        MatrixXf id1= id*(0.00001);
        MatrixXf ncxx = cxx+id1;
        MatrixXf ncyy = cyy+id1;


        //Taking inverse.
        MatrixXf icxx= ncxx.inverse();
        MatrixXf icyy= ncyy.inverse();
  cout <<endl<< "cxx "<<endl<<cxx<<endl;
  cout <<endl<< "cxy"<<endl<<cxy<<endl;
  cout <<endl<< "cyx"<<endl<<cyx<<endl;
  cout <<endl<< "cyy"<<endl<<cyy<<endl;
  


        cout<<endl<<endl<<"icxx"<< endl<<icxx<<endl;
        cout<<"icyy"<<endl<<icyy<<endl;


        MatrixXf finpro = (icxx*cxy*icyy*cyx);
        cout<<endl<<"finpro"<<endl<<finpro<<endl;
        //Eigenvalue decomposition
   EigenSolver<MatrixXf> eig(finpro);
   //Printing the eigenvalues.
   int r=0;
   cout<<"eigen values are  "<<endl;
   for(r=0;r<nrowsa;r++)
   {
           complex<double> eigenval1 = sqrt(eig.eigenvalues()[r]);
                   cout<<eigenval1<<endl;
   }


    //Printing the eigenvectors.
        cout << "eigenvec"<<endl<< eig.eigenvectors() << endl<<endl<<endl;


   


}
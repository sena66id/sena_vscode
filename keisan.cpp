#include<iostream>
#include<fstream>
#include<vector>
#include<string>
#include<sstream>
#include<cmath>
#include<algorithm>
#include <Eigen/Core>
#include <Eigen/LU>

using namespace Eigen;
using namespace std;
void export_vtu(const std::string &file, vector<vector<double>> node, vector<vector<int>> element, vector<double> C)
{
    FILE *fp;
  if ((fp = fopen(file.c_str(), "w")) == NULL)
  {
    cout << file << " open error" << endl;
    exit(1);
  }
  fprintf(fp, "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt32\">\n");
  fprintf(fp, "<UnstructuredGrid>\n");
  fprintf(fp, "<Piece NumberOfPoints= \"%d\" NumberOfCells= \"%d\" >\n", node.size(), element.size());
  fprintf(fp, "<Points>\n");
  int offset = 0;
  fprintf(fp, "<DataArray type=\"Float64\" Name=\"Position\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%d\"/>\n",offset);
  offset += sizeof(int) + sizeof(double) * node.size() * 3;
  fprintf(fp, "</Points>\n");
  fprintf(fp, "<Cells>\n");
  fprintf(fp, "<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n");
  for (int i = 0; i < element.size(); i++){
    for (int j = 0; j < element[i].size(); j++) fprintf(fp, "%d ", element[i][j]);
    fprintf(fp, "\n");
  }
  fprintf(fp, "</DataArray>\n");
  fprintf(fp, "<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n");
  int num = 0;
  for (int i = 0; i < element.size(); i++)
  {
    num += element[i].size();
    fprintf(fp, "%d\n", num);
  }
  fprintf(fp, "</DataArray>\n");
  fprintf(fp, "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
  for (int i = 0; i < element.size(); i++) fprintf(fp, "%d\n", 5);
  fprintf(fp, "</DataArray>\n");
  fprintf(fp, "</Cells>\n");
  fprintf(fp, "<PointData>\n");
  fprintf(fp, "<DataArray type=\"Float64\" Name=\"pressure[Pa]\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%d\"/>\n",offset);
  offset += sizeof(int) + sizeof(double) * node.size();
  fprintf(fp, "</PointData>\n");
  fprintf(fp, "<CellData>\n");
  fprintf(fp, "</CellData>\n");
  fprintf(fp, "</Piece>\n");
  fprintf(fp, "</UnstructuredGrid>\n");
  fprintf(fp, "<AppendedData encoding=\"raw\">\n");
  fprintf(fp, "_");
  fclose(fp);
  fstream ofs;
  ofs.open(file.c_str(), ios::out | ios::app | ios_base::binary);
  double *data_d = new double[node.size()*3];
  num = 0;
  int size=0;
  for (int ic = 0; ic < node.size(); ic++){
    data_d[num] = node[ic][0];
    num++;
    data_d[num] = node[ic][1];
    num++;
    data_d[num] = 0.0;
    num++;
  }
  size=sizeof(double)*node.size()*3;
  ofs.write((char *)&size, sizeof(size));
  ofs.write((char *)data_d, size);
  num=0;
  for (int ic = 0; ic < node.size(); ic++){
      data_d[num]   = C[ic];
      num++;
  }
  size=sizeof(double)*node.size();
  ofs.write((char *)&size, sizeof(size));
  ofs.write((char *)data_d, size);
  delete data_d;
  ofs.close();
  if ((fp = fopen(file.c_str(), "a")) == NULL)
  {
    cout << file << " open error" << endl;
    exit(1);
  }
  fprintf(fp, "\n</AppendedData>\n");
  fprintf(fp, "</VTKFile>\n");
  fclose(fp);
}

MatrixXd K(1301,1301);
//double K[1301][1301];
VectorXd L(1301);
VectorXd U(1301);
int main()
{
    string str;
    ifstream ifs("node.dat");
    //vector<double> t(2) は double t[2]と同じ
    //vector<vector<double>> t(2, vector<double>(2))は double t[2][2]と同じ
    vector<vector<double>> x;
    while(getline(ifs,str)){
        istringstream ss(str);
        string tmp;
        vector<double> tmp_x;
        for(int j=0; j<3; j++){
            getline(ss, tmp, ' ');
            tmp_x.push_back(stod(tmp));
        }
        x.push_back(tmp_x);
    }
    ifs.close();
    ifs.open("element.dat");
    vector<vector<int>> element;
    while(getline(ifs,str)){
        istringstream ss(str);
        string tmp;
        vector<int> tmp_element;
        for(int j=0; j<4; j++){
            getline(ss, tmp, ' ');
            if(j==0) continue;
            tmp_element.push_back(stoi(tmp));
        }
        element.push_back(tmp_element);
    }
    ifs.close();
    /*ifs.close();
    vector<double> C(x.size(),0.0);
    double minimum = 100000.0;
    for(int i=0; i<x.size(); i++){
        minimum=min(minimum,x[i][0]);
    }
    ofstream ofs("boundary.dat");
    for(int i=0; i<x.size(); i++){
        if(fabs(x[i][0]-minimum)<0.000001){
            ofs << i << endl;
            C[i] = 1.0;
        }
    }
    ofs.close();*/
    
    vector<double> C1(x.size(),0.0);
    double minimum = 100000.0;
    for(int i=0; i<x.size(); i++){
        minimum=min(minimum,x[i][0]);
    }
    ofstream ofs1("boundary_left.dat");
    for(int i=0; i<x.size(); i++){
        if(fabs(x[i][0]-minimum)<0.000001){
            ofs1 << i << endl;
            C1[i] = 1.0;
        }
    }
    ofs1.close();

    vector<double> C2(x.size(),0.0);
    double maximum = 0.001;
    for(int i=0; i<x.size(); i++){
        maximum=max(maximum,x[i][0]);
    }
    ofstream ofs2("boundary_right.dat");
    for(int i=0; i<x.size(); i++){
        if(fabs(x[i][0]-maximum)<0.000001){
            ofs2 << i << endl;
            C2[i] = 1.0;
        }
    }
    ofs2.close();
    //export_vtu("test.vtu", x, element, C1);

    ifs.open("boundary_left.dat");
    vector<int> boundary_left;
    while(getline(ifs,str)){
        istringstream ss(str);
        string tmp;
        getline(ss, tmp, ' ');
        boundary_left.push_back(stoi(tmp));
    }
    ifs.close();

    ifs.open("boundary_right.dat");
    vector<int> boundary_right;
    while(getline(ifs,str)){
        istringstream ss(str);
        string tmp;
        getline(ss, tmp, ' ');
        boundary_right.push_back(stoi(tmp));
    }
    ifs.close();


       
    




    //double K[x.size()][x.size()]={0};

    //for(int i=0;i<x.size();i++){
        //for(int j=0;j<x.size();j++){
          //  K[i][j]=0.0;

      //  }
    

  //  }

    
    
    double dxdr[2][2];
    double drdx[2][2];
    double dNdr[2][3];
    double dNdx[2][3];
    double Ke[3][3]; double w=1;
    dNdr[0][0]=-1;
    dNdr[0][1]=1;
    dNdr[0][2]=0;
    dNdr[1][0]=-1;
    dNdr[1][1]=0;
    dNdr[1][2]=1;
    double dxdr1[2],dxdr2[2],dxdr3[2];
    //cout<<x[element[0][0]][0]<<endl;
    
    
    

    for(int ic=0;ic<1301;ic++ ){
        L(ic)=0.0;
        U(ic)=0.0;
        for(int in=0;in<1301;in++){
            //K[in][ic]=0.0;
            K(in,ic)=0.0;
        }

    }


    for(int ic=0;ic<1301;ic++){
        for(int j=0;j<3;j++){
            for(int i=0; i<2;i++){

                dxdr1[i]=0.0;
                dxdr2[i]=0.0;
                
                dNdx[i][j]=0.0;


            }
            for(int n=0; n<3;n++){
                Ke[j][n]=0.0;

            }
        }
        

            for(int j=0;j<2;j++){
                for(int l=0;l<3;l++){
                    dxdr1[j] +=dNdr[j][l]*x[element[ic][l]][0];
                    dxdr2[j] +=dNdr[j][l]*x[element[ic][l]][1];
                }
                //dxdr1[j] +=dNdr(j,0)*x[element[0][j]][0];
                //dxdr2[j] +=dNdr(j,1)*x[element[0][j]][1];    
                /*cout<<dxdr1[j]<<endl;
                cout<<dxdr2[j] << endl;*/
            }        
            
        
        double J=dxdr1[0]*dxdr2[1]-dxdr1[1]*dxdr2[0];
        //cout<<J<<endl;
        drdx[0][0]=(dxdr2[1])/J;drdx[0][1]=-(dxdr1[1])/J;
        drdx[1][0]=-(dxdr2[0])/J;drdx[1][1]=(dxdr1[0])/J;
        
        for(int i=0;i<2;i++){
            for(int j=0;j<2;j++){
                //cout<<drdx[i][j];
                //cout<<" ";
            }
        //cout<<""
        }
        //cout<<"\n";
        
        for(int i=0;i<3 ;i++){
            for(int j=0; j<2 ;j++){

                for(int p=0;p<2;p++){
                    dNdx[j][i] +=dNdr[p][i]*drdx[p][j];
                    
                }
                /*cout<<dNdx[j][i];
                cout <<"  ";*/
            }
            //cout<<endl;
        
        }
        for(int i=0;i<3;i++){
            for(int j=0;j<3;j++){
        
                for(int p=0;p<2;p++ ){
            
                    Ke[i][j]+=dNdx[p][i]*dNdx[p][j]*J;
                
                }
                /*cout<<Ke[i][j];
                cout<<" ";*/
                
            }

            //cout<<"\n"<<endl;

        }

        /*cout<<"\n"<<endl;
        cout<<"\n"<<endl;*/
        for(int i=0; i<3;i++){
            for(int j=0;j<3;j++ ){

                K(element[ic][i],element[ic][j])+=Ke[i][j];
            }
            
        }
        //for(int in=0;in<1301;in++){
            //cout<<K();
           // cout<<" ";
          
        //}
      //cout<<"\n";
        



        /*MatrixXd a = MatrixXd::Ones(2,2);

        a(0,0) = 2;

        std::cout << a << std::endl;*/

        //for(int i=0;i<x.size();i++)

        

    }

    

    
            
        
        ofstream outputfile1("K.txt");
        outputfile1<< K ;
        outputfile1.close();

    for(int i=0;i<boundary_right.size();i++){
        for(int j=0;j<1301;j++){
            K(boundary_right[i],j)=0.0;
        }
    K(boundary_right[i],boundary_right[i])=1.0;
    L(boundary_right[i])=1.0;

    }    

    for(int i=0;i<boundary_left.size();i++){
        for(int j=0;j<1301;j++){
            K(boundary_left[i],j)=0.0;
        }
    K(boundary_left[i],boundary_left[i])=1.0;
    L(boundary_left[i])=100.0;

    }    
    

    ofstream outputfile2("L.txt");
    outputfile2<<L;
    outputfile2.close();
    
    
    
    FullPivLU<MatrixXd> lu(K);
    U = lu.solve(L);
    ofstream outputfile3("U.txt");
    outputfile3 << U;
    outputfile3.close();
    
    vector<double> p(x.size());
    for(int i=0; i<x.size();i++){
        p[i]=U(i);
    }

   





        return 0;
}





    





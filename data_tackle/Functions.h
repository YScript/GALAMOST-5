/*
GALAMOST - GPU-Accelerated Large-Scale Molecular Simulation Toolkit
Version 5
COPYRIGHT
	GALAMOST Copyright (c) (2020) G.U.G.D.
LICENSE
	This program is a free software: you can redistribute it and/or 
	modify it under the terms of the GNU General Public License. 
	This program is distributed in the hope that it will be useful, 
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANT ABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
	See the General Public License v3 for more details.
	You should have received a copy of the GNU General Public License
	along with this program. If not, see <http://www.gnu.org/licenses/>.
DISCLAIMER
	The authors of GALAMOST do not guarantee that this program and its 
	derivatives are free from error. In no event shall the copyright 
	holder or contributors be liable for any indirect, incidental, 
	special, exemplary, or consequential loss or damage that results 
	from its use. We also have no responsibility for providing the 
	service of functional extension of this program to general users.
USER OBLIGATION 
	If any results obtained with GALAMOST are published in the scientific 
	literature, the users have an obligation to distribute this program 
	and acknowledge our efforts by citing the paper "Y.-L. Zhu et al.,
	J. Comput. Chem. 2013,34, 2197-2211" in their article.
CORRESPONDENCE
	Dr. You-Liang Zhu, 
	Email: ylzhu@galamost.com
*/

#include <cuda_runtime.h>
#include "MolInfo.h"

#ifndef __FUNCTIONS_H__
#define __FUNCTIONS_H__

class Function
	{
    public:
        Function() {};
        virtual ~Function() {};	
		void add(mst_reader* build, MolInfo* mol)
			{
			m_build = build;
			m_mol = mol;
			};
		virtual void compute(){};
		mst_reader* m_build;		
		MolInfo* m_mol;
	};

//--- case 1
class Rg2 : public Function
    {
    public:
        Rg2(std::string filename): Function()
			{
			m_Nf =0;
			m_file.open(filename.c_str(), ios::app);
			if (!m_file.good())
				{
				cerr << endl << "***Error! Error opening dump file " << filename << endl << endl;
				throw runtime_error("Error Rg2 dump");
				}
			}
		virtual ~Rg2() 
			{
			for(unsigned int i=0; i<m_av_rg2.size();i++)
				cout<<"Mol"<<i<<" rg2= "<<m_av_rg2[i]/double(m_Nf)<<endl;
			};
		virtual void compute();		
    private:
		std::ofstream m_file;
		std::vector<double> m_av_rg2;
		unsigned int m_Nf;
	};
//--- case 2	
class Ed2 : public Function
    {
    public:
        Ed2(std::string filename): Function()
			{
			m_Nf =0;			
			m_file.open(filename.c_str(), ios::app);
			if (!m_file.good())
				{
				cerr << endl << "***Error! Error opening dump file " << filename << endl << endl;
				throw runtime_error("Error Ed2 dump");
				}
			}
		virtual ~Ed2()
			{
			for(unsigned int i=0; i<m_av_ed2.size();i++)
				cout<<"Mol"<<i<<" ed2 = "<<m_av_ed2[i]/double(m_Nf)<<endl;
			};		

		virtual void compute();			
    private:
		std::ofstream m_file;
		std::vector<double> m_av_ed2;
		unsigned int m_Nf;		
	};
//--- case 3
class RDF : public Function
    {
    public:
        RDF(std::string filename): Function()
			{
			m_file.open(filename.c_str(), ios_base::out);
			if (!m_file.good())
				{
				cerr << endl << "***Error! Error opening dump file " << filename << endl << endl;
				throw runtime_error("Error RDF dump");
				}
			m_maxbin =100;
			m_block_size = 256;
			m_gpu_id = 0;
			m_rdf.resize(m_maxbin);
			m_r.resize(m_maxbin);
			m_Nf=0;
			m_rmax=0.0;
			m_exclusion_mol = false;
			m_exclusion_list = false;
			m_exclusion_angle = false;
			m_exclusion_bond = false;			
			}
		virtual ~RDF() 
			{
			for(unsigned int i=0; i<m_rdf.size(); i++)
				m_file<<m_r[i]<<"  "<<m_rdf[i]/double(m_Nf)<<"\n";
			m_file.close();
			}
		void setMaxbin(unsigned int maxbin)
			{
			m_maxbin=maxbin;
			m_rdf.clear();
			m_rdf.resize(m_maxbin);
			m_r.clear();
			m_r.resize(m_maxbin);
			}
		void setRmax(double rmax)
			{
			m_rmax = rmax;
			}
		void setGPU(unsigned int gpu_id)
			{
			m_gpu_id = gpu_id;
			cudaSetDevice(m_gpu_id);
			}
		void setBondEx(bool bond_ex)
			{
			m_exclusion_bond = bond_ex;
			}			
		void setAngleEx(bool angle_ex)
			{
			m_exclusion_angle = angle_ex;
			}	
		void setMolEx(bool mol_ex)
			{
			m_exclusion_mol = mol_ex;
			}				
		virtual void compute();			
    private:
		std::ofstream m_file;
		unsigned int m_maxbin;
		unsigned int m_block_size;
		int m_gpu_id;
		std::vector<double> m_rdf;
		std::vector<double> m_r;		
		unsigned int m_Nf;
		double m_rmax;
		bool m_exclusion_mol;
		bool m_exclusion_list;
		bool m_exclusion_angle;
		bool m_exclusion_bond;
		unsigned int* d_n_exclusion;
		unsigned int* d_exclusion_list;
		unsigned int* d_mol_id_per_particle;
	};
//--- case 4
class Bond_distr : public Function
    {
    public:
        Bond_distr(std::string filename): Function()
			{
			m_file.open(filename.c_str(), ios_base::out);
			if (!m_file.good())
				{
				cerr << endl << "***Error! Error opening dump file " << filename << endl << endl;
				throw runtime_error("Error Bond_distr dump");
				}
			m_Npot = 2001;
			m_Nf=0;			
			}
		virtual ~Bond_distr() 
			{
			for(unsigned int i=0; i<m_Nb;i++)
				{
				m_file<<m_bondMap[i]<<endl;
				for(unsigned int j=0;j<m_Npot;j++)
					{
					if(m_distr[i*m_Npot+j]>0)
						{
						double value = m_distr[i*m_Npot+j]/double(m_Nf);
						m_file<<double(j)*m_delt<<"  "<<value/m_delt<<"\n";
						}
					}
				cout<<"The averaged length of bond "<<m_bondMap[i]<<" is "<<m_bond_lenth[i]/double(m_Nf)<<endl;					
				}
			m_file.close();				
			}
		void setNpot(unsigned int npot)
			{
			m_Npot=npot;
			}
		virtual void compute();			
    private:
		std::ofstream m_file;
		double m_rcut;
		unsigned int m_Npot;
		unsigned int m_Nf;
		unsigned int m_Nb;
		double m_delt;		
		std::vector<double> m_distr;
		std::vector<double> m_bond_lenth;		
		std::vector<std::string> m_bondMap;
	};
//--- case 5	
class Angle_distr : public Function
    {
    public:
        Angle_distr(std::string filename): Function()
			{
			m_file.open(filename.c_str(), ios_base::out);
			if (!m_file.good())
				{
				cerr << endl << "***Error! Error opening dump file " << filename << endl << endl;
				throw runtime_error("Error Angle_distr dump");
				}
			m_Npot = 2001;	
			m_Nf=0;				
			}
		virtual ~Angle_distr() 
			{
			for(unsigned int i=0; i<m_Nb;i++)
				{
				m_file<<m_angleMap[i]<<endl;
				for(unsigned int j=0;j<m_Npot;j++)
					{
					if(m_distr[i*m_Npot+j]>0)
						{
						double value = m_distr[i*m_Npot+j]/double(m_Nf);
						m_file<<double(j)*m_delt<<"  "<<value/m_delt<<"\n";
						}
					}
				cout<<"The averaged radian of angle "<<m_angleMap[i]<<" is "<<m_angle_radian[i]/double(m_Nf)<<endl;							
				}
			m_file.close();						
			}
		void setNpot(unsigned int npot)
			{
			m_Npot=npot;
			}
		virtual void compute();			
    private:
		std::ofstream m_file;
		unsigned int m_Npot;
		unsigned int m_Nf;
		unsigned int m_Nb;
		double m_delt;
		std::vector<double > m_distr;
		std::vector<double> m_angle_radian;		
		std::vector<std::string> m_angleMap;		
	};
//--- case 6
class Dihedral_distr : public Function
    {
    public:
        Dihedral_distr(std::string filename): Function()
			{
			m_file.open(filename.c_str(), ios_base::out);
			if (!m_file.good())
				{
				cerr << endl << "***Error! Error opening dump file " << filename << endl << endl;
				throw runtime_error("Error Dihedral_distr dump");
				}
			m_Npot = 2001;	
			m_Nf=0;					
			}
		virtual ~Dihedral_distr() 
			{
			for(unsigned int i=0; i<m_Nb;i++)
				{
				m_file<<m_dihedralMap[i]<<endl;
				for(unsigned int j=0;j<m_Npot;j++)
					{
					if(m_distr[i*m_Npot+j]>0)
						{
						double value = m_distr[i*m_Npot+j]/double(m_Nf);
						m_file<<double(j)*m_delt<<"  "<<value/m_delt<<"\n";
						}
					}
				cout<<"The averaged radian of dihedral "<<m_dihedralMap[i]<<" is "<<m_dihedral_radian[i]/double(m_Nf)<<endl;					
				}
			m_file.close();						
			}
		void setNpot(unsigned int npot)
			{
			m_Npot=npot;
			}
		virtual void compute();			
    private:		
		std::ofstream m_file;
		unsigned int m_Npot;
		unsigned int m_Nf;
		unsigned int m_Nb;
		double m_delt;	
		std::vector<double > m_distr;
		std::vector<double> m_dihedral_radian;			
		std::vector<std::string> m_dihedralMap;
	};
//--- case 7
class StressTensor : public Function
    {
    public:
        StressTensor(std::string filename): Function()
			{
			m_file.open(filename.c_str(), ios::app);
			if (!m_file.good())
				{
				cerr << endl << "***Error! Error opening dump file " << filename << endl << endl;
				throw runtime_error("Error StressTensor dump");
				}
			m_Nf=0;
			m_bondex=true;
			m_bodyex=true;
			m_diameter = true;
			m_rcut=3.0;
			m_delt=0.0;
			}
		virtual ~StressTensor() {};
		void setParam();	
		virtual void compute();
		void setBondEx(bool bondex)
			{
			m_bondex = bondex;
			}
		void setBodyEx(bool bodyex)
			{
			m_bodyex = bodyex;
			}
		void setDiameterConsider(bool diameter)
			{
			m_diameter = diameter;
			}
    private:
		std::ofstream m_file;
		unsigned int m_Nf;
		std::vector<vec> m_pparams;
		std::vector<vec> m_bparams;
		double m_rcut;
		double m_delt;
		bool m_bondex;
		bool m_bodyex;
		bool m_diameter;
	};
//--- case 8
class Density : public Function
    {
    public:
        Density(std::string filename): Function()
			{
			m_file.open(filename.c_str(), ios::app);
			if (!m_file.good())
				{
				cerr << endl << "***Error! Error opening dump file " << filename << endl << endl;
				throw runtime_error("Error Density dump");
				}
			m_Nf=0;
			}
		virtual ~Density() {};

		virtual void compute();			
    private:
		std::ofstream m_file;
		unsigned int m_Nf;
	};	
//--- case 9
class Reimage : public Function
    {
    public:
        Reimage(): Function()
			{
			m_solute_center = false;
			m_label_free_particle = false;
			m_unwrap_molecule = true;
			m_keep_molecule_center_in_box = false;
			m_remove_bond_cross_box = false;
			m_remove_image = false;
			m_body_keep = false;
			m_target_type = "";
			m_shiftx = 0.0;
			m_shifty = 0.0;
			m_shiftz = 0.0;
			m_nprecision = 10;
			m_nhead = 7;
			m_add_image_to_pos = true;
			m_Nf=0;
			m_sp = ios::beg;
			m_file="xxxxxxxx";
			}
		virtual ~Reimage() {};
		void setShiftX(double shiftx)
			{
			m_shiftx = shiftx;
			}
		void setShiftY(double shifty)
			{
			m_shifty = shifty;
			}
		void setShiftZ(double shiftz)
			{
			m_shiftz = shiftz;
			}
		void setLabelFreeParticle(std::string target_type)
			{
			m_target_type = target_type;
			m_label_free_particle = true;
			}
		void setUnwrapMolecule(bool unwrap_molecule)
			{
			m_unwrap_molecule = unwrap_molecule;
			}
		void setMoleculeCenterInBox(bool keep_molecule_center_in_box)
			{
			m_keep_molecule_center_in_box = keep_molecule_center_in_box;
			}	
		void setRemoveImage(bool remove_image)
			{
			m_remove_image = remove_image;
			m_add_image_to_pos=false;
			}
		void setRemoveBondCrossBox(bool remove_bond_cross_box)
			{
			m_remove_bond_cross_box = remove_bond_cross_box;
			}
		void addImageToPos(bool add_image_to_pos)
			{
			m_add_image_to_pos = add_image_to_pos;
			}
		void setBodyKeep(bool body_keep)
			{
			m_body_keep = body_keep;
			}				
		virtual void compute();			
    private:
		bool m_solute_center;
		bool m_label_free_particle;
		bool m_unwrap_molecule;
		bool m_keep_molecule_center_in_box;
		bool m_remove_image;
		bool m_remove_bond_cross_box;
		bool m_add_image_to_pos;
		bool m_body_keep;
		double m_shiftx;
		double m_shifty;
		double m_shiftz;
		unsigned int m_nprecision;
		unsigned int m_nhead;
		std::string m_target_type;
		unsigned int m_Nf;
		ifstream::pos_type m_sp;
		std::string m_file;
	};	
//--- case 10
class MSD : public Function
    {
    public:
        MSD(std::string filename): Function()
			{
			m_file.open(filename.c_str(), ios_base::out);
			if (!m_file.good())
				{
				cerr << endl << "***Error! Error opening dump file " << filename << endl << endl;
				throw runtime_error("Error MSD dump");
				}
			m_Nf=0;
			}
		virtual ~MSD() {};
		void setParam();
		void setDirection(std::string direction)
			{
			m_direction=direction;
			}
		virtual void compute();			
    private:
		unsigned int m_Nf;
		std::string m_direction;
		std::ofstream m_file;
		std::vector<vec> m_pos_offset;
		std::vector<vec> m_pos_cm;		
	};
//--- case 11
class RDFCM : public Function
    {
    public:
        RDFCM(std::string filename): Function()
			{
			m_file.open(filename.c_str(), ios_base::out);
			if (!m_file.good())
				{
				cerr << endl << "***Error! Error opening dump file " << filename << endl << endl;
				throw runtime_error("Error RDFCM dump");
				}
			m_maxbin =100;
			m_block_size = 256;
			m_gpu_id = 0;
			m_rdf.resize(m_maxbin);
			m_r.resize(m_maxbin);			
			m_Nf=0;
			m_rmax = 0.0;
			m_exclusion_mol = false;
			m_exclusion_list = false;
			m_exclusion_angle = false;
			m_exclusion_bond = false;		
			}
		virtual ~RDFCM() 
			{
			for(unsigned int i=0; i<m_rdf.size(); i++)
				m_file<<m_r[i]<<"  "<<m_rdf[i]/double(m_Nf)<<"\n";
			m_file.close();
			}
		void setMaxbin(unsigned int maxbin)
			{
			m_maxbin=maxbin;
			m_rdf.clear();
			m_rdf.resize(m_maxbin);
			m_r.clear();
			m_r.resize(m_maxbin);
			}
		void setRmax(double rmax)
			{
			m_rmax = rmax;
			}			
		void setGPU(unsigned int gpu_id)
			{
			m_gpu_id = gpu_id;
			cudaSetDevice(m_gpu_id);
			}
		virtual void compute();			
    private:
		std::ofstream m_file;
		unsigned int m_maxbin;
		unsigned int m_block_size;
		int m_gpu_id;
		std::vector<double> m_rdf;
		std::vector<double> m_r;		
		unsigned int m_Nf;
		double m_rmax;
		bool m_exclusion_mol;
		bool m_exclusion_list;
		bool m_exclusion_angle;
		bool m_exclusion_bond;	
		unsigned int* d_n_exclusion;
		unsigned int* d_exclusion_list;
		unsigned int* d_mol_id_per_particle;		
	};	
//--- case 12
class MSDCM : public Function
    {
    public:
        MSDCM(std::string filename): Function()
			{
			m_file.open(filename.c_str(), ios_base::out);
			if (!m_file.good())
				{
				cerr << endl << "***Error! Error opening dump file " << filename << endl << endl;
				throw runtime_error("Error MSDCM dump");
				}
			m_Nf=0;
			}
		virtual ~MSDCM();
		virtual void compute();	
		void setDirection(std::string direction)
			{
			m_direction=direction;
			}		
    private:
		unsigned int m_Nf;
		std::string m_direction;		
		std::ofstream m_file;
		std::vector<std::vector<vec> > m_pos_all;
		std::vector<unsigned int > m_mol_type_id;
		unsigned int m_n_kind_mol;	
	};

//--- case 13
class Entanglement : public Function
    {
    public:
        Entanglement(std::string filename): Function()
			{
			m_file.open(filename.c_str(), ios_base::out);
			if (!m_file.good())
				{
				cerr << endl << "***Error! Error opening dump file " << filename << endl << endl;
				throw runtime_error("Error Entanglement dump");
				}
			m_Nf=0;
			m_rcut_max = 0.0;
			m_totalnum = 0;
			m_fdelt = 1.0;
			m_fmax = 10000.0;
			m_fdistr.resize((unsigned int)(m_fmax/m_fdelt));
			}
		virtual ~Entanglement() 
			{
			for(unsigned int i=0; i<m_fdistr.size();i++)
				{
				if(m_fdistr[i]>0)
					{
					double distr = double(m_fdistr[i])/(double(m_totalnum)*m_fdelt);
					m_file<<(double(i)+0.5)*m_fdelt<<"  "<<distr<<"  "<<m_fdistr[i]<<endl;
					}
				}
			m_file.close();
			};
		void setParam();			
		virtual void compute();			
    private:
		unsigned int m_Nf;
		std::ofstream m_file;
		std::vector<double> m_rcutsq;
		std::vector<unsigned int> m_fdistr;
		unsigned int m_totalnum;
		double m_fdelt;
		double m_fmax;
		std::vector<vec> m_params;		
		double m_rcut_max;
	};	
//--- case 14
class STRFAC : public Function
    {
    public:
        STRFAC(std::string filename): Function()
			{
			m_file.open(filename.c_str(), ios_base::out);
			if (!m_file.good())
				{
				cerr << endl << "***Error! Error opening dump file " << filename << endl << endl;
				throw runtime_error("Error STRFAC dump");
				}
			m_Kmax=80;
			m_LKmax=0;
			m_MKmax=0;
			m_NKmax=0;			
			m_Nf=0;
			m_Ntypes=0;
			m_block_size = 256;
			m_gpu_id = 0;
			m_L = 0.0;
			m_Qmax = 0.0;
			m_direction = "XYZ";
			m_2D=false;
			}
		virtual ~STRFAC() 
			{
			unsigned int Ksqmax = m_Kmax*m_Kmax;
			if(m_2D)
				{
				m_file<<"qx"<<"  "<<"qy"<<"  ";
				for(unsigned int typi =0; typi <m_Ntypes; typi++)
					{
					for(unsigned int typj =typi; typj <m_Ntypes; typj++)
						{
						m_file<<m_type_map[typi]+"-"+m_type_map[typj]+"_r"<<"  "<<m_type_map[typi]+"-"+m_type_map[typj]+"_i"<<"  ";
						}
					}
				m_file<<"\n";							

				for(unsigned int ii =1; ii<(unsigned int)(m_LKmax+1); ii++)
					{
					unsigned int i = m_LKmax + 1 - ii;
					double qi =double(i)*M_PI*2.0/m_L;						
					for(unsigned int jj =0; jj<(unsigned int)(2*m_MKmax+1); jj++)
						{
						unsigned int j = 2*m_MKmax -jj;
						double qj =double(int(j)-int(m_MKmax))*M_PI*2.0/m_L;
						m_file<<-qi<<"   "<<-qj<<"   ";
						unsigned int icn =0;
						for(unsigned int typi =0; typi < m_Ntypes; typi++)
							{	
							for(unsigned int typj =typi; typj < m_Ntypes; typj++)
								{	
								unsigned int idi = i+icn*(m_LKmax+1);
								unsigned int idj = j+icn*(2*m_MKmax+1);
								m_file<<m_sqreal2D[idi][idj]/double(m_Nf)<<"   "<<m_sqimage2D[idi][idj]/double(m_Nf)<<"  ";
								icn +=1;
								}
							}
						m_file<<"\n";	
						}
					m_file<<"\n";						
					}

				for(unsigned int i =0; i<(unsigned int)(m_LKmax+1); i++)
					{
					double qi =double(i)*M_PI*2.0/m_L;						
					for(unsigned int j =0; j<(unsigned int)(2*m_MKmax+1); j++)
						{
						double qj =double(int(j)-int(m_MKmax))*M_PI*2.0/m_L;
						m_file<<qi<<"   "<<qj<<"   ";
						unsigned int icn =0;
						for(unsigned int typi =0; typi < m_Ntypes; typi++)
							{	
							for(unsigned int typj =typi; typj < m_Ntypes; typj++)
								{	
								unsigned int idi = i+icn*(m_LKmax+1);
								unsigned int idj = j+icn*(2*m_MKmax+1);
								if(i==0&&int(j)<m_MKmax)
									idj = 2*m_MKmax-j+icn*(2*m_MKmax+1);
								m_file<<m_sqreal2D[idi][idj]/double(m_Nf)<<"   "<<m_sqimage2D[idi][idj]/double(m_Nf)<<"  ";
								icn +=1;
								}
							}
						m_file<<"\n";	
						}
					m_file<<"\n";	
					}						
				}
			else
				{
				m_file<<"q"<<"  ";
				for(unsigned int typi =0; typi <m_Ntypes; typi++)
					{
					for(unsigned int typj =typi; typj <m_Ntypes; typj++)
						{
						m_file<<m_type_map[typi]+"-"+m_type_map[typj]+"_r"<<"  "<<m_type_map[typi]+"-"+m_type_map[typj]+"_i"<<"  ";
						}
					}
				m_file<<"\n";					
					
				for(unsigned int j =0; j<Ksqmax; j++)
					{
					unsigned int icn = 0;
					float sumx=0.0;
					float sumy=0.0;
					for(unsigned int typi =0; typi < m_Ntypes; typi++)
						{
						for(unsigned int typj =typi; typj < m_Ntypes; typj++)
							{								
							unsigned int id = icn*Ksqmax + j;
							sumx += m_sqreal[id]/double(m_Nf);
							sumy += m_sqimage[id]/double(m_Nf);
							icn +=1;
							}
						}
					if(sumx!=0.0||sumy!=0.0)
						{						
						unsigned int icn = 0;
						double q =sqrt(double(j+1))*M_PI*2.0/m_L;
						m_file<<q<<"  ";
						for(unsigned int typi =0; typi < m_Ntypes; typi++)
							{
							for(unsigned int typj =typi; typj < m_Ntypes; typj++)
								{								
								unsigned int id = icn*Ksqmax + j;
								m_file<<m_sqreal[id]/double(m_Nf)<<"  "<<m_sqimage[id]/double(m_Nf)<<"  ";
								icn +=1;
								}
							}
						m_file<<"\n";
						}					
					}						
				}
			m_file.close();
			}
		void setQmax(float qmax)
			{
			m_Qmax=qmax;
			}
		void setGPU(unsigned int gpu_id)
			{
			m_gpu_id = gpu_id;
			cudaSetDevice(m_gpu_id);
			}
		void setDeltaq(float deltaq)
			{
			m_L = M_PI*2.0/deltaq;
			}
		void setDirection(string direction)
			{
			m_direction=direction;
			}	
		void set2D(bool d2)
			{
			m_2D=d2;
			}				
		virtual void compute();			
    private:
		std::ofstream m_file;
		int m_Kmax;
		int m_LKmax;
		int m_MKmax;
		int m_NKmax;		
		float m_Qmax;
		std::vector<double> m_sqreal;
		std::vector<double> m_sqimage;
		std::vector<std::vector<double> > m_sqreal2D;
		std::vector<std::vector<double> > m_sqimage2D;

		std::vector<string> m_type_map;
		unsigned int m_Ntypes;
		unsigned int m_Nf;
		unsigned int m_gpu_id;
		unsigned int m_block_size;
		float m_L;
		bool m_2D;
		string m_direction;		
	};
//--- case 15
class DOMAINSIZE : public Function
    {
    public:
        DOMAINSIZE(std::string filename): Function()
			{
			m_file.open(filename.c_str(), ios_base::out);
			if (!m_file.good())
				{
				cerr << endl << "***Error! Error opening dump file " << filename << endl << endl;
				throw runtime_error("Error DOMAINSIZE dump");
				}
			m_Kmax=20;
			m_Nf=0;
			m_block_size = 256;
			m_gpu_id = 0;
			m_qc = 0.4690;     //scf=0.3708;0.4690 //dpd=0.6188;0.5809
			}
		virtual ~DOMAINSIZE() 
			{
			m_file.close();
			}
		void setKmax(unsigned int kmax)
			{
			m_Kmax=kmax;
			}
		void setQc(float qc)
			{
			m_qc=qc;
			}
		
		void setGPU(unsigned int gpu_id)
			{
			m_gpu_id = gpu_id;
			cudaSetDevice(m_gpu_id);
			}
		virtual void compute();			
    private:
		std::ofstream m_file;
		unsigned int m_Kmax;
		unsigned int m_Nf;
		unsigned int m_gpu_id;
		unsigned int m_block_size;
		float m_qc;
	};
//--- case 16
class DSTRFAC : public Function
    {
    public:
        DSTRFAC(std::string filename): Function()
			{
			m_file.open(filename.c_str(), ios_base::out);
			if (!m_file.good())
				{
				cerr << endl << "***Error! Error opening dump file " << filename << endl << endl;
				throw runtime_error("Error DSTRFAC dump");
				}
			m_Kmax=vec_int(0, 0, 0);
			m_q=7.0;
			m_Nf=0;
			}
		virtual ~DSTRFAC() {};
		void setParam();
		virtual void compute();	
		void setKmax(unsigned int kmax)
			{
			m_Kmax=vec_int(kmax, kmax, kmax);
			}
		void setQ(double q)
			{
			m_q=q;
			}
    private:
		unsigned int m_Nf;
		vec_int m_Kmax;
		double m_q;
		std::vector<vec_int> m_qvec;	
		std::ofstream m_file;
		std::vector<vec> m_pos_offset;
		std::vector<vec> m_pos_cm;		
	};
//--- case 17
class ConfigCheck : public Function
    {
    public:
        ConfigCheck(std::string filename): Function()
			{
			m_file.open(filename.c_str(), ios_base::out);
			if (!m_file.good())
				{
				cerr << endl << "***Error! Error opening dump file " << filename << endl << endl;
				throw runtime_error("Error ConfigCheck dump");
				}
			m_Nf=0;
			m_bondex=true;
			m_angleex=true;
			m_dihedralex=true;
			m_bodyex=true;
			m_rcut=2.0;
			}
		virtual ~ConfigCheck() {};
		virtual void compute();
		void setBondEx(bool bondex)
			{
			m_bondex = bondex;
			}
		void setAngleEx(bool angleex)
			{
			m_angleex = angleex;
			}
		void setDihedralEx(bool dihedralex)
			{
			m_dihedralex = dihedralex;
			}
		void setBodyEx(bool bodyex)
			{
			m_bodyex = bodyex;
			}
		void setRcut(double rcut)
			{
			m_rcut = rcut;
			}
    private:
		unsigned int m_Nf;	
		std::ofstream m_file;
		double m_rcut;
		bool m_bondex;
		bool m_bodyex;
		bool m_angleex;
		bool m_dihedralex;
	};
	
//--- case 18
class RDFBetweenTypes : public Function
    {
    public:
        RDFBetweenTypes(std::string filename): Function()
			{
			m_file.open(filename.c_str(), ios_base::out);
			if (!m_file.good())
				{
				cerr << endl << "***Error! Error opening dump file " << filename << endl << endl;
				throw runtime_error("Error RDFBetweenTypes dump");
				}
			m_maxbin =100;
			m_block_size = 256;
			m_gpu_id = 0;
			m_rdf.resize(m_maxbin);
			m_r.resize(m_maxbin);
			m_Nf=0;
			m_rmax = 0.0;
			m_exclusion_mol = false;
			m_exclusion_list = false;
			m_exclusion_angle = false;
			m_exclusion_bond = false;	
			}
		virtual ~RDFBetweenTypes() 
			{
			m_file<<"r"<<"  ";
			for(unsigned int typi =0; typi <m_Ntype; typi++)
				{
				for(unsigned int typj =typi; typj <m_Ntype; typj++)
					{
					m_file<<m_type_map[typi]+"-"+m_type_map[typj]<<"  ";
					}
				}
			m_file<<"\n";
			for (unsigned int bin = 0; bin < m_maxbin; bin++ )
				{ 
				double r = m_r[bin];
				m_file<<r<<"  ";
				for(unsigned int typi =0; typi <m_Ntype; typi++)
					{
					for(unsigned int typj =typi; typj <m_Ntype; typj++)
						{
						double g = m_rdf[(typi*m_Ntype+typj)*m_maxbin+bin]/double(m_Nf);
						m_file<<g<<"  ";
						}  
					}
				m_file<<"\n";
				}
			m_file.close();
			}
		void setMaxbin(unsigned int maxbin)
			{
			m_maxbin=maxbin;
			m_rdf.clear();
			m_rdf.resize(m_maxbin);
			m_r.clear();
			m_r.resize(m_maxbin);
			}
		void setRmax(double rmax)
			{
			m_rmax = rmax;
			}				
		void setGPU(unsigned int gpu_id)
			{
			m_gpu_id = gpu_id;
			cudaSetDevice(m_gpu_id);
			}
		void setBondEx(bool bond_ex)
			{
			m_exclusion_bond = bond_ex;
			}			
		void setAngleEx(bool angle_ex)
			{
			m_exclusion_angle = angle_ex;
			}	
		void setMolEx(bool mol_ex)
			{
			m_exclusion_mol = mol_ex;
			}			
		virtual void compute();			
    private:
		std::ofstream m_file;
		unsigned int m_maxbin;
		unsigned int m_block_size;
		int m_gpu_id;
		std::vector<double> m_rdf;
		std::vector<double> m_r;		
		unsigned int m_Nf;
		std::vector<std::string> m_type_map;
		unsigned int m_Ntype;
		double m_rmax;
		bool m_exclusion_mol;
		bool m_exclusion_list;
		bool m_exclusion_angle;
		bool m_exclusion_bond;
		unsigned int* d_n_exclusion;
		unsigned int* d_exclusion_list;
		unsigned int* d_mol_id_per_particle;		
	};

class MSTfileConversion : public Function
    {
    public:
        MSTfileConversion(): Function()
			{
			m_lammps=false;
			m_gromacs = false;
			m_xml=false;
			m_Nf = 0;
			m_nprecision=10;
			m_nhead = 7;
			}
		virtual ~MSTfileConversion() 
			{
			}
		void setLammps(bool lammps)
			{
			m_lammps = lammps;
			}
		void setGromacs(bool gromacs)
			{
			m_gromacs = gromacs;
			}
		void setXML(bool xml)
			{
			m_xml = xml;
			}		
		virtual void compute();			
    private:
		std::ofstream m_file;
		bool m_lammps;
		bool m_gromacs;
		bool m_xml;
		unsigned int m_Nf;
		unsigned int m_nprecision;
		unsigned int m_nhead;		
	};	
	
class PatchToParticle : public Function
    {
    public:
        PatchToParticle(): Function()
			{
			m_nprecision=10;
			m_nhead = 7;
			m_separatedis = 0.1;
			m_scale = 0.97;
			m_filter_sphere = false;
			m_Nf = 0;
			}
		virtual ~PatchToParticle() 
			{
			}
		void setSeparateDis(double sd)
			{
			m_separatedis = sd;
			}
		void setFilterSphere(bool fs)
			{
			m_filter_sphere = fs;
			}
		void setPatchParticleScale(double sc)
			{
			m_scale = sc;
			}			
		virtual void compute();			
    private:
		std::ofstream m_file;
		unsigned int m_nprecision;
		unsigned int m_nhead;
		double m_separatedis;
		double m_scale;
		bool m_filter_sphere;
		unsigned int m_Nf;
	};		
	
#endif


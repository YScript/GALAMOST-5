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

#include "mst_reader.h"

std::vector<std::string> split(std::string str, std::string pattern)
	{
	istringstream i_stream(str);
	std::vector<std::string> q;
	std::string s;
	while (i_stream>>s) {q.push_back(s);}
    return q;
	}

mst_reader::mst_reader()
    {
    m_timestep = 0;
    m_dimension = 3;
	m_box = BoxSize(0.0, 0.0, 0.0);
	m_num_particles = 0;
	m_mst_read = false;
	m_invariant_data = false;
	m_variant_data = false;
	m_if_trajectory = false;
	m_read_indicator = { {"bond", false}, {"angle", false}, {"box", false}, {"position", false}, {"type", false}, 
	                     {"image", false}, {"mass", false}, {"velocity", false}, {"force", false}, {"virial", false} };

	m_sp = ios::beg;
    m_fname = "XXXXXXXXX";
	m_object_name = "mst_reader";
    }

void mst_reader::reset_params()
	{
	m_num_particles_read = false;
	m_timestep_read = false;
	m_dimension_read = false;
	m_position_read = false;
	m_type_read = false;
	m_image_read = false;
	m_mass_read = false;
	m_velocity_read = false;
	m_bond_read = false;
	m_angle_read = false;
	m_box_read = false;
	}
	
void mst_reader::clear_data()
	{
	if (!m_read_indicator["position"])				
		m_pos.clear();
	if (!m_read_indicator["type"])		
		m_type.clear();
	if (!m_read_indicator["image"])		
		m_image.clear();
	if (!m_read_indicator["mass"])			
		m_mass.clear();
	if (!m_read_indicator["velocity"])			
		m_vel.clear();
	if (!m_read_indicator["force"])			
		m_force.clear();
	if (!m_read_indicator["virial"])			
		m_virial.clear();
	if (!m_read_indicator["bond"])
		m_bonds.clear();
	if (!m_read_indicator["angle"])		
		m_angles.clear();
	if (!m_read_indicator["box"])				
		m_box=BoxSize(0.0, 0.0, 0.0);
	}


unsigned int mst_reader::getNDimensions() const
    {
    return (unsigned int)m_dimension;
    }

unsigned int mst_reader::getNParticles() const
    {
    return m_num_particles;
    }

unsigned int mst_reader::getNParticleTypes() const
    {
    return (unsigned int)m_type_map.size();
    }

BoxSize mst_reader::getBox() const
    {
    return m_box;
    }

unsigned int mst_reader::getTimeStep() const
    {
    return m_timestep;
    }

bool mst_reader::readDataFromMST(const string &fname)
    {
	m_if_changed_np = false;
	reset_params();
	bool read_frame = false;
	ifstream file;
	file.open(fname.c_str());	

	if (!file.good())
		{
		cerr << endl << "Unable to open file " << fname.c_str() << endl << endl;
		throw runtime_error("Error reading mst file");
		}

	if (m_fname == fname)
		{
		file.seekg(m_sp);
		}
	else
		{
		// cout<<"read '"<<fname.c_str()<<"'"<<endl;
		m_fname = fname;		
		m_mst_read = false;
		m_invariant_data = false;
		m_variant_data = false;	
		m_if_trajectory = false;
		for (std::map<std::string, bool>::iterator it = m_read_indicator.begin(); it != m_read_indicator.end();it++)
			it->second = false;
		
		file.seekg(0, ios::beg);
		clear_data();		
		}
		
	std::string line;

	while(getline(file, line))
		{
		if (line.find("mst_version") != line.npos && line.find("1.0") != line.npos )
			{
			// cout<<"read mst file with version 1.0"<<endl;
			m_mst_read = true;
			continue;
			}			
		// cout<<line<<endl;
		if (m_mst_read)
			{
			if (line.find("invariant_data") != line.npos)
				{
				m_invariant_data = true;
				m_variant_data = false;						
				continue;
				}
				
			if (line.find("variant_data") != line.npos)
				{
				m_invariant_data = false;
				m_variant_data = true;	
				continue;
				}
				
			if (line.find("frame") != line.npos)
				{
				// check node 
					{
					istringstream parser;
					parser.str(line);
					std::string name;
					int frame_id = -1;
					parser >> name;
					parser >> frame_id;
					
					if (name != "frame"||frame_id==-1)
						throw runtime_error("Error! mst file with wrong format for the indicator of 'frame'");
					
					m_if_trajectory = true;
					}

				if (read_frame)
					{
					file.seekg(-line.size()-1, ios::cur);
					m_sp=file.tellg();	
					break;
					}
					
				if (m_variant_data)
					clear_data();
				else
					throw runtime_error("Error! mst files with multiple frames without the label of 'variant_data'");
				read_frame = true;
				continue;
				}				

			if (line.find("num_particles") != line.npos)
				{
				reset_params();
				m_num_particles_read=true;
				continue;
				}
				
			if (line.find("timestep") != line.npos)
				{
				reset_params();
				m_timestep_read=true;
				continue;
				}

			if (line.find("dimension") != line.npos)
				{
				reset_params();
				m_dimension_read=true;
				continue;
				}				
		
			if (line.find("position") != line.npos)
				{		
				reset_params();			
				m_position_read=true;
				if (m_invariant_data)
					m_read_indicator["position"] = true;					
				continue;
				}

			if (line.find("type") != line.npos)
				{
				reset_params();			
				m_type_read=true;
				if (m_invariant_data)
					m_read_indicator["type"] = true;					
				continue;
				}

			if (line.find("image") != line.npos)
				{				
				reset_params();
				m_image_read=true;
				if (m_invariant_data)
					m_read_indicator["image"] = true;
				continue;
				}

			if (line.find("mass") != line.npos)
				{
				reset_params();
				m_mass_read=true;
				if (m_invariant_data)
					m_read_indicator["mass"] = true;						
				continue;
				}
				
			if (line.find("velocity") != line.npos)
				{
				reset_params();
				m_velocity_read=true;
				if (m_invariant_data)
					m_read_indicator["velocity"] = true;
				continue;
				}
				
			if (line.find("bond") != line.npos)
				{
				reset_params();
				m_bond_read=true;	
				if (m_invariant_data)
					m_read_indicator["bond"] = true;
				continue;
				}

			if (line.find("angle") != line.npos)
				{				
				reset_params();
				m_angle_read=true;
				if (m_invariant_data)
					m_read_indicator["angle"] = true;
				continue;
				}

			if (line.find("box") != line.npos)
				{
				reset_params();
				m_box_read=true;
				if (m_invariant_data)
					m_read_indicator["box"] = true;
				continue;
				}
				
			// read data				
			std::vector<std::string> line_array = split(line, "	");
			istringstream parser;
			parser.str(line);

			if (m_num_particles_read && line_array.size() >= 1)
				{
				unsigned int np;
				parser >> np;
				if (np!=m_num_particles)
					m_if_changed_np =  true;			
				m_num_particles = np;

				}
				
			if (m_timestep_read && line_array.size() >= 1)
				parser >> m_timestep;
				
			if (m_dimension_read && line_array.size() >= 1)
				parser >> m_dimension;		
				
			if (m_box_read && line_array.size() >= 3)
				{
				double lx, ly, lz;
				parser>> lx >> ly >> lz;
				m_box = BoxSize(lx, ly, lz);
				}
				
			if (m_position_read && line_array.size() >= 3)
				{
				double px, py, pz;				
				parser>> px >> py >> pz;
				// cout<<px<<" "<<py<<" "<<pz<<endl;				
				m_pos.push_back(vec(px, py, pz));
				}

			if (m_type_read && line_array.size() >= 1)
				{
				// cout<<line_array[0]<<endl;
				string type;
				parser>> type;
				m_type.push_back(getTypeId(type));
				}
				
			if (m_image_read && line_array.size() >= 3)
				{
				int ix, iy, iz;
				parser>> ix >> iy >> iz;
				m_image.push_back(vec_int(ix, iy, iz));
				}
				
			if (m_mass_read && line_array.size() >= 1)
				{
				double mass;
				parser>> mass;
				m_mass.push_back(mass);
				}
				
			if (m_velocity_read && line_array.size() >= 3)
				{
				double vx, vy, vz;
				parser>> vx >> vy >> vz;
				m_vel.push_back(vec(vx, vy, vz));
				}
				
			if (m_bond_read && line_array.size() >= 3)
				{
				string name;
				unsigned int a;
				unsigned int b;
				parser>> name >> a >> b;
				m_bonds.push_back(Bond(name, a, b, getBondTypeId(name)));
				}
		
			if (m_angle_read && line_array.size() >= 4)
				{
				string name;
				unsigned int a;
				unsigned int b;
				unsigned int c;
				parser>> name >> a >> b >> c;
				m_angles.push_back(Angle(name, a, b, c, getAngleTypeId(name)));
				}		
			}
		}
		
	if(!m_mst_read)
		{
        cerr << endl
             << "***Error! This is not a MST file"
             << endl << endl;
        throw runtime_error("Error extracting data from galamost_mst file");			
		}		
		
    if (m_box.lx==0.0&&m_box.ly==0.0&&m_box.lz==0.0)
        {
        cerr << endl
             << "***Error! A box is required to define simulation box"
             << endl << endl;
        throw runtime_error("Error extracting data from galamost_mst file");
        }
		

	if(m_dimension==2&&m_box.lz>0.0)
		{
        cerr << "***Error! two dimensions of the simulation box should be with Lz = 0.0, the Lz = " << m_box.lz <<" in galamost_mst files"<< endl << endl;
        throw runtime_error("Error extracting data from galamost_mst file");
		}
	
	if(m_dimension==3&&m_box.lz<1.0e-6)
		{
        cerr << "***Error! Lz = 0.0 should be with two dimensions of the simulation, three dimensions defind in galamost_mst files "<< endl << endl;
        throw runtime_error("Error extracting data from galamost_mst file");
		}
	
    if (m_pos.size() == 0)
        {
        cerr << endl << "***Error! No particles defined in position node" << endl << endl;
        throw runtime_error("Error extracting data from galamost_mst file");
        }
		
    if (m_type.size() == 0)
        {
        cerr << endl << "***Error! No particles defined in type node" << endl << endl;
        throw runtime_error("Error extracting data from galamost_mst file");
        }
		
		
    if (m_pos.size() != m_num_particles)
        {
        cerr << endl << "***Error! " << m_pos.size() << " positions != " << m_num_particles
             << " the number of particles" << endl << endl;
        throw runtime_error("Error extracting data from galamost_mst file");
        }
		
    if (m_type.size() != m_num_particles)
        {
        cerr << endl << "***Error! " << m_type.size() << " types != " << m_num_particles
             << " the number of particles" << endl << endl;
        throw runtime_error("Error extracting data from galamost_mst file");
        }	
		
    if (m_molecule.size() != 0 && m_molecule.size() != m_num_particles)
        {
        cerr << endl << "***Error! " << m_molecule.size() << " molecule != " << m_num_particles
             << " the number of particles" << endl << endl;
        throw runtime_error("Error extracting data from galamost_mst file");
        }
		
    if (m_vel.size() != 0 && m_vel.size() != m_num_particles)
        {
        cerr << endl << "***Error! " << m_vel.size() << " velocities != " << m_num_particles
             << " the number of particles" << endl << endl;
        throw runtime_error("Error extracting data from galamost_mst file");
        }
		
    if (m_mass.size() != 0 && m_mass.size() != m_num_particles)
        {
        cerr << endl << "***Error! " << m_mass.size() << " masses != " << m_num_particles
             << " the number of particles" << endl << endl;
        throw runtime_error("Error extracting data from galamost_mst file");
        }
		
    if (m_diameter.size() != 0 && m_diameter.size() != m_num_particles)
        {
        cerr << endl << "***Error! " << m_diameter.size() << " diameters != " << m_num_particles
             << " the number of particles" << endl << endl;
        throw runtime_error("Error extracting data from galamost_mst file");
        }
		
    if (m_image.size() != 0 && m_image.size() != m_num_particles)
        {
        cerr << endl << "***Error! " << m_image.size() << " images != " << m_num_particles
             << " the number of particles" << endl << endl;
        throw runtime_error("Error extracting data from galamost_mst file");
        }
		
    if ( m_type.size() != 0 && m_type.size() != m_num_particles)
        {
        cerr << endl << "***Error! " << m_type.size() << " type values != " << m_num_particles
             << " the number of particles" << endl << endl;
        throw runtime_error("Error extracting data from galamost_mst file");
        }
		
    if (m_body.size() != 0 && m_body.size() != m_num_particles)
        {
        cerr << endl << "***Error! " << m_body.size() << " body values != " << m_num_particles
             << " the number of particles" << endl << endl;
        throw runtime_error("Error extracting data from galamost_mst file");
        }	
		
    if (m_charge.size() != 0 && m_charge.size() != m_num_particles)
        {
        cerr << endl << "***Error! " << m_charge.size() << " charge values != " << m_num_particles
             << " the number of particles" << endl << endl;
        throw runtime_error("Error extracting data from galamost_mst file");
        }
		
	//-check bonds, angles and dihedrals
	for (unsigned int i=0; i<m_bonds.size();i++)
		{
		Bond bi = m_bonds[i];
		if(bi.a>=m_num_particles||bi.b>=m_num_particles)
			{
			cerr << endl << "***Error! bond '" << bi.type <<" "<<bi.a<<" "<<bi.b<< "' with the particle index not in the range [0, N-1], with N = "<< m_num_particles<<" ! " << endl << endl;
			throw runtime_error("Error extracting data from galamost_mst file");				
			}
		}
		
	for (unsigned int i=0; i<m_angles.size();i++)
		{
		Angle ai = m_angles[i];
		if(ai.a>=m_num_particles||ai.b>=m_num_particles||ai.c>=m_num_particles)
			{
			cerr << endl << "***Error! angle '" << ai.type <<" "<<ai.a<<" "<<ai.b<<" "<<ai.c<< "' with the particle index not in the range [0, N-1], with N = "<< m_num_particles<<" ! " << endl << endl;
			throw runtime_error("Error extracting data from galamost_mst file");				
			}
		}	
		
	for (unsigned int i=0; i<m_dihedrals.size();i++)
		{
		Dihedral di = m_dihedrals[i];
		if(di.a>=m_num_particles||di.b>=m_num_particles||di.c>=m_num_particles||di.d>=m_num_particles)
			{
			cerr << endl << "***Error! dihedral '" << di.type <<" "<<di.a<<" "<<di.b<<" "<<di.c<<" "<<di.d<< "' with the particle index not in the range [0, N-1], with N = "<< m_num_particles<<" ! " << endl << endl;
			throw runtime_error("Error extracting data from galamost_mst file");				
			}
		}

	for (unsigned int i=0; i<m_vsites.size();i++)
		{
		Dihedral di = m_vsites[i];
		if(di.a>=m_num_particles||di.b>=m_num_particles||di.c>=m_num_particles||di.d>=m_num_particles)
			{
			cerr << endl << "***Error! vsite '" << di.type <<" "<<di.a<<" "<<di.b<<" "<<di.c<<" "<<di.d<< "' with the particle index not in the range [0, N-1], with N = "<< m_num_particles<<" ! " << endl << endl;
			throw runtime_error("Error extracting data from galamost_mst file");				
			}
		}

	if (file.eof())
		return false;
	
	return read_frame;
    }
	
void mst_reader::outPutInfo()
	{
    // notify the user of what we have accomplished
    cout <<"--- galamost mst file read summary" << endl;
    cout <<" "<< getNParticles() << " particles at timestep " << m_timestep << endl;
    if (m_image.size() > 0)
        cout <<" "<< m_image.size() << " images" << endl;
    if (m_vel.size() > 0)
        cout <<" "<< m_vel.size() << " velocities" << endl;
    if (m_mass.size() > 0)
        cout <<" "<< m_mass.size() << " masses" << endl;
    if (m_diameter.size() > 0)
        cout <<" "<< m_diameter.size() << " diameters" << endl;
    cout <<" "<< getNParticleTypes() <<  " particle types" << endl;
    if (m_body.size() > 0)
        cout <<" "<< m_body.size() << " particle body values" << endl; 	
    if (m_bonds.size() > 0)
        cout <<" "<< m_bonds.size() << " bonds" << endl;
    if (m_angles.size() > 0)
        cout <<" "<< m_angles.size() << " angles" << endl;
    if (m_dihedrals.size() > 0)
        cout <<" "<< m_dihedrals.size() << " dihedrals" << endl;
    if (m_charge.size() > 0)
        cout <<" "<< m_charge.size() << " charges" << endl;
    if (m_orientation.size() > 0)
        cout <<" "<< m_orientation.size() << " orientations" << endl;
    if (m_quaternion.size() > 0)
        cout <<" "<< m_quaternion.size() << " quaternions" << endl;		
    if (m_molecule.size() > 0)
        cout <<" "<< m_molecule.size() << " molecules" << endl;	
	if(m_mass.size()==0)
		{
		m_mass.resize(m_pos.size());
		for(unsigned int i=0; i<m_mass.size();i++)
			m_mass[i] = 1.0;
        cout <<" "<<" set mass to be 1.0 by default!" << endl;			
		}		
	}


unsigned int mst_reader::getTypeId(const std::string& name) 
    {
    for (unsigned int i = 0; i < m_type_map.size(); i++)
        {
        if (m_type_map[i] == name)
            return i;
        }
    m_type_map.push_back(name);
    return (unsigned int)m_type_map.size()-1;
    }

unsigned int mst_reader::getBondTypeId(const std::string& name)
    {
    for (unsigned int i = 0; i < m_bond_type_map.size(); i++)
        {
        if (m_bond_type_map[i] == name)
            return i;
        }
    m_bond_type_map.push_back(name);
    return (unsigned int)m_bond_type_map.size()-1;
    }


unsigned int mst_reader::getAngleTypeId(const std::string& name)
    {
    for (unsigned int i = 0; i < m_angle_type_map.size(); i++)
        {
        if (m_angle_type_map[i] == name)
            return i;
        }
    m_angle_type_map.push_back(name);
    return (unsigned int)m_angle_type_map.size()-1;
    }


unsigned int mst_reader::getDihedralTypeId(const std::string& name)
    {
    for (unsigned int i = 0; i < m_dihedral_type_map.size(); i++)
        {
        if (m_dihedral_type_map[i] == name)
            return i;
        }
    m_dihedral_type_map.push_back(name);
    return (unsigned int)m_dihedral_type_map.size()-1;
    }

unsigned int mst_reader::getVsiteTypeId(const std::string& name)
    {
    for (unsigned int i = 0; i < m_vsite_type_map.size(); i++)
        {
        if (m_vsite_type_map[i] == name)
            return i;
        }
    m_vsite_type_map.push_back(name);
    return (unsigned int)m_vsite_type_map.size()-1;
    }

unsigned int mst_reader::getNBondTypes()  const
    {
    return (unsigned int)m_bond_type_map.size();
    }

unsigned int mst_reader::getNAngleTypes() const 
    {
    return (unsigned int)m_angle_type_map.size();
    }

unsigned int mst_reader::getNDihedralTypes() const
    {
    return (unsigned int)m_dihedral_type_map.size();
    }

unsigned int mst_reader::getNVsiteTypes() const
    {
    return (unsigned int)m_vsite_type_map.size();
    }


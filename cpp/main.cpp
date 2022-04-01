#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <iostream>
#include <fstream>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/convert_nef_polyhedron_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_slicer.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef CGAL::Nef_polyhedron_3<Kernel> Nef_polyhedron;
typedef CGAL::Surface_mesh<Kernel::Point_3> Surface_mesh;

template <class HDS>
struct Polyhedron_builder : public CGAL::Modifier_base<HDS> {
    std::vector<Point> vertices; // type: Kernel::Point_3, for EACH SHELL
    std::vector<std::vector<unsigned long>> faces; // INDEX for vertices in EACH SHELL

    Polyhedron_builder() {}
    void operator()(HDS& hds) {
        CGAL::Polyhedron_incremental_builder_3<HDS> builder(hds, true);
        std::cout << "constructing polyhedron -- ";
        std::cout << "building surface with " << vertices.size() << " vertices and " << faces.size() << " faces" << '\n';

        builder.begin_surface(vertices.size(), faces.size());
        for (auto const& vertex : vertices) builder.add_vertex(vertex);
        for (auto const& face : faces) builder.add_facet(face.begin(), face.end());
        builder.end_surface();
    }
};


/*
 * this function is for split the face with delimiter "\\"
 */
std::vector<std::string> split(const std::string& str, const std::string& delim) {
    std::vector<std::string> res;
    if(str.empty()) return res;
    char * strs = new char[str.length() + 1] ;
    strcpy(strs, str.c_str());
    char * d = new char[delim.length() + 1];
    strcpy(d, delim.c_str());
    char *p = strtok(strs, d);
    while(p) {
        std::string s = p;
        res.push_back(s);
        p = strtok(NULL, d);
    }
    return res;
}

/*
 * this function is to write test output file, as obj format
 */
static void output_OBJ_file(const std::vector<Point>& vertices, const std::vector<std::vector<unsigned long>> faces,std::string& fname) {

    std::ofstream output_file;

    output_file.open(fname, std::ios::out | std::ios::trunc);
    output_file << std::fixed;

    // print data into the file
    for(const auto& ver:vertices)
        output_file<<"v " << ver.x()<<" "<< ver.y()<<";" << ver.z()<<std::endl;
    for(const auto &face:faces)
    {
        output_file<<"f "<<face[0]<<" "<<face[1]<<" "<<face[2]<<std::endl;
    }
    // close file
    output_file.close();
}

/*
 * this function is to read ifc file and save polyhedron as a vector
 */
std::vector<Polyhedron> read_ifc(std::string& fname) {

    std::vector<Polyhedron> polyhedrons;
    std::string INTER_PATH = "../../data/"; // the data folder, ie D: .../data
    std::string filename = INTER_PATH + fname;
    std::cout << "-- loading obj shell: " << filename << '\n';

    // read obj
    std::string line;
    std::ifstream file(filename);
    if (!file.is_open()) { std::cerr << "file open failed! " << '\n'; }

    Polyhedron_builder<Polyhedron::HalfedgeDS> polyhedron_builder; // construct polyhedron_builder

    std::vector<Point> all_vertices;
    std::vector<std::vector<unsigned long>> face_v_indices;
    unsigned long  vertex_mark = 0;
    std::vector<Point> cur_vertices;
    // process each line in the obj file

    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string word;
        iss >> word;
        if (line.empty())
        {
            continue; // skip empty lines:
        }
        if (word == "v") {
            std::vector<double> coordinates;
            while (iss >> word) coordinates.push_back(std::stof(word));
            if (coordinates.size() == 3) {
                cur_vertices.emplace_back(coordinates[0], coordinates[1], coordinates[2]);
                all_vertices.emplace_back(coordinates[0], coordinates[1], coordinates[2]);
            }
            else cur_vertices.emplace_back();
        }
        else if (word == "f") {
            std::vector<unsigned long> current_face;
            while (iss >> word) {
                std::string str = split(word, "//")[0];
                unsigned long num = stoi(str);
                current_face.emplace_back(num);
            }
            face_v_indices.emplace_back(current_face);
        }
        else if(word == "haha")
        {
            std::cout<<"hahahahaha"<<std::endl;
        }

        else if (word == "s") {
            if(cur_vertices.empty()) continue;

            else if (!cur_vertices.empty() && !face_v_indices.empty()) {

                std::vector<Point> unique_vertices;
                // this function is to remove redundant vertices
                for(unsigned long i=0 ; i<cur_vertices.size() ; ++i)
                {
                    std::cout<<"current vertices"<<cur_vertices[i].x()<<", "<<cur_vertices[i].y()<<", "<<cur_vertices[i].z()<<std::endl;
                    if(unique_vertices.empty()) unique_vertices.emplace_back(cur_vertices[i]);
                    else
                    {
                        bool mark= false;
                        for(unsigned long j=0;j<unique_vertices.size();++j)
                        {
                            if(i != j)
                            {
                                if(cur_vertices[i].x() == unique_vertices[j].x()&&cur_vertices[i].y() == unique_vertices[j].y() && cur_vertices[i].z() == unique_vertices[j].z() )
                                {
                                    mark= true;
                                }
                            }
                        }
                        if(!mark) {
                            unique_vertices.emplace_back(cur_vertices[i]);
                        }
                    }
                }

                for(auto & face:face_v_indices)
                {
                    for(auto & vertex: face)
                    {

                        for(unsigned long i = 0; i < unique_vertices.size(); ++i){
                            if(unique_vertices[i].x()==all_vertices[vertex-1].x() && unique_vertices[i].y()==all_vertices[vertex-1].y() && unique_vertices[i].z()==all_vertices[vertex-1].z())
                            {
                                vertex=i;
                                break;
                            }
                        }
                    }
                }
                for (auto &coordinates: unique_vertices) {
                    polyhedron_builder.vertices.emplace_back(
                            Point(
                                    coordinates[0], // x
                                    coordinates[1], // y
                                    coordinates[2]) // z
                    );
                }
                for (auto &face: face_v_indices) {
                    polyhedron_builder.faces.emplace_back(face);
                }
                Polyhedron polyhedron;
                polyhedron.delegate(polyhedron_builder);
                polyhedrons.emplace_back(polyhedron);

                std::cout << "build polyhedron now" << '\n';

                std::cout << "done" << '\n';
                face_v_indices.clear();
                cur_vertices.clear();
                polyhedron_builder.vertices.clear();
                polyhedron_builder.faces.clear();
            }
        }
    }
    if(!cur_vertices.empty())
    {
        std::vector<Point> unique_vertices;
        // this function is to remove redundant vertices
        for(unsigned long i=0 ; i<cur_vertices.size() ; ++i)
        {
            std::cout<<"current vertices"<<cur_vertices[i].x()<<", "<<cur_vertices[i].y()<<", "<<cur_vertices[i].z()<<std::endl;
            if(unique_vertices.empty()) unique_vertices.emplace_back(cur_vertices[i]);
            else
            {
                bool mark= false;
                for(unsigned long j=0;j<unique_vertices.size();++j)
                {
                    if(i != j)
                    {
                        if(cur_vertices[i].x() == unique_vertices[j].x()&&cur_vertices[i].y() == unique_vertices[j].y() && cur_vertices[i].z() == unique_vertices[j].z() )
                        {
                            mark= true;
                            std::cout<<"i is "<<i<<", j is "<<j<<std::endl;
                        }
                    }
                }
                if(!mark) {
                    unique_vertices.emplace_back(cur_vertices[i]);
                    std::cout<<"add to unique vertices: "<<cur_vertices[i].x()<<", "<<cur_vertices[i].y()<<", "<<cur_vertices[i].z()<<std::endl;
                    std::cout<<"unique vertices size: "<<unique_vertices.size()<<std::endl;
                }
            }
        }

        for(auto & face:face_v_indices)
        {
            std::cout<<"current face is "<<face[0]<<", "<<face[1]<<", "<<face[2]<<std::endl;
            for(auto & vertex: face)
            {

                for(unsigned long i = 0; i < unique_vertices.size(); ++i){
                    if(unique_vertices[i].x()==all_vertices[vertex-1].x() && unique_vertices[i].y()==all_vertices[vertex-1].y() && unique_vertices[i].z()==all_vertices[vertex-1].z())
                    {
                        std::cout<<"i is "<<i<<", vertex is "<<vertex<<std::endl;
                        vertex=i;
                        break;
                    }
                }
            }
        }
        for (auto &coordinates: unique_vertices) {
            polyhedron_builder.vertices.emplace_back(
                    Point(
                            coordinates[0], // x
                            coordinates[1], // y
                            coordinates[2]) // z
            );
        }
        for (auto &face: face_v_indices) {
            polyhedron_builder.faces.emplace_back(face);
        }
        Polyhedron polyhedron;
        polyhedron.delegate(polyhedron_builder);
        polyhedrons.emplace_back(polyhedron);

        std::cout << "build polyhedron now" << '\n';

        std::cout << "done" << '\n';
        face_v_indices.clear();
        cur_vertices.clear();
        polyhedron_builder.vertices.clear();
        polyhedron_builder.faces.clear();
    }
    return polyhedrons;
}

/*
 * construct polyhedron for each obj file(each obj file stands for one shell)
 */
Polyhedron construct_polyhedron_each_shell_obj(std::string &fname) {

        std::string INTER_PATH = "../../data"; // the data folder, ie D: .../data
        std::string filename = INTER_PATH + fname;
        std::cout << "-- loading obj shell: " << filename << '\n';

        // read obj
        std::string line;
        std::ifstream file(filename);
        if (!file.is_open()) { std::cerr << "file open failed! " << '\n'; }

        std::vector<double> coordinates; // store xyz coordinates of each vertex line
        std::vector<unsigned long> face_v_indices; // store face-vertex indices in each face line

        Polyhedron_builder<Polyhedron::HalfedgeDS> polyhedron_builder; // construct polyhedron_builder

        // process each line in the obj file
        while (std::getline(file, line)) {
            if (line.empty()) {
                continue; // skip empty lines:
            }

            // use markers below to mark the type of each line
            bool v_flag(false); // entitled "v"
            bool f_flag(false); // entitled "f"

            // help to process elements in each line
            std::stringstream ss(line);
            std::string field; // element in each line
            std::string::size_type sz; // NOT std::size_t, size of std::string

            // for each element(field, type: string) in one line
            while (std::getline(ss, field, ' ')) {

                // get the type of current line
                if (field == "f") {
                    f_flag = true;
                    continue;
                } else if (field == "v") {
                    v_flag = true;
                    continue;
                }

                // process the current line
                if (v_flag) {
                    // process xyz coordinates
                    coordinates.emplace_back(std::stod(field, &sz));
                } else if (f_flag) {
                    face_v_indices.emplace_back(std::stoi(field, &sz) - 1);
                }

            } // end while: process each element in one line

            // process each vertex (if it's a vertex line)
            if (!coordinates.empty() && coordinates.size() == 3) {
                polyhedron_builder.vertices.emplace_back(
                        Point(
                                coordinates[0], // x
                                coordinates[1], // y
                                coordinates[2]) // z
                );
            }
            coordinates.clear();

            // add face-v indices
            if (!face_v_indices.empty()) {
                polyhedron_builder.faces.emplace_back(face_v_indices);
            }
            face_v_indices.clear();

        } // end while: each line in the file

        // construct polyhedron ------------------------------------------------------------------
        Polyhedron polyhedron;
        polyhedron.delegate(polyhedron_builder);
        std::cout << "done" << '\n';
        return polyhedron;
}

/*
 * this function is to create big nef polyhedron
 */
Nef_polyhedron create_nef(std::vector<Polyhedron> polyhedron_list){
    Nef_polyhedron nef;
    for(auto &polyhedra:polyhedron_list){
        Nef_polyhedron nef_po=Nef_polyhedron(polyhedra);
        nef=nef+nef_po;
    }
    return nef;
}

/*
 * this function is to write test output file, as OFF format
 */
void output_OFF_file(const Nef_polyhedron& nef_p,const std::string& f) {
    Surface_mesh output;
    CGAL::convert_nef_polyhedron_to_polygon_mesh(nef_p, output);
    std::ofstream out;
    //std::string f="out_mesh.OFF";
    out.open(f);
    out << output;
    out.close();
}

/*
 * this is shell explorer
 */
struct Shell_explorer {
    std::vector<Point> vertices;
    std::vector<std::vector<unsigned long>> faces;

    void visit(Nef_polyhedron::Vertex_const_handle v) {}
    void visit(Nef_polyhedron::Halfedge_const_handle he) {}
    void visit(Nef_polyhedron::SHalfedge_const_handle she) {}
    void visit(Nef_polyhedron::SHalfloop_const_handle shl) {}
    void visit(Nef_polyhedron::SFace_const_handle sf) {}

    void visit(Nef_polyhedron::Halffacet_const_handle hf) {

    }
};


int main()
{
    std::string fname = "test2.obj";
    std::vector<Polyhedron> polyhedrons_list = read_ifc(fname);
    Nef_polyhedron big_nef;
    for(auto & polyhedron:polyhedrons_list)
    {
        Nef_polyhedron temp_nef=Nef_polyhedron (polyhedron);
        big_nef+=temp_nef;
    }
    output_OFF_file(big_nef,"output_off.OFF");

    Nef_polyhedron::Volume_const_iterator current_volume;
    CGAL_forall_volumes(current_volume, big_nef) {
        Nef_polyhedron::Shell_entry_const_iterator current_shell;
        CGAL_forall_shells_of(current_shell, current_volume) {
            Shell_explorer se;
            Nef_polyhedron::SFace_const_handle sface_in_shell(current_shell);
            big_nef.visit_shell_objects(sface_in_shell, se);

        }
    }

    return 0;
}

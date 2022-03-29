#include <iostream>
#include <fstream>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Nef_polyhedron_3.h>


typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef CGAL::Nef_polyhedron_3<Kernel> Nef_polyhedron;

struct Face {
    std::vector<int> vertices; // indices in a vector of points
};

class Shell {
public:
    std::vector<Face> faces;
};

class Object {
public:
    std::string id;
    std::vector<Shell> shells;
};

// Next, load each shell into a CGAL Polyhedron_3 using the Polyhedron_incremental_builder_3.
// need to construct class to use CGAL
template <class HDS>

struct Polyhedron_builder : public CGAL::Modifier_base<HDS> {
    std::vector<Point> vertices;
    std::vector<std::vector<Point>> faces;

    Polyhedron_builder() = default;
    void operator()(HDS& hds) {
        CGAL::Polyhedron_incremental_builder_3<HDS> builder(hds, true);
        std::cout << "building surface with " << vertices.size() << " vertices and " << faces.size() << " faces" << std::endl;

        builder.begin_surface(vertices.size(), faces.size());
        for (auto const &vertex: vertices) builder.add_vertex(vertex);
        for (auto const &face: faces) builder.add_facet(face.begin(), face.end());
        builder.end_surface();
    }
};

bool vertex_in_list(std::vector<float> current_vertex,std::vector<Point> vertices) {
    for(auto& vertex: vertices){
        if(vertex[0]==current_vertex[0] && vertex[1]==current_vertex[1] && vertex[2] == current_vertex[2]) {
            return true;
        }
    }
    return false;
}

// this function is for split the face with delimiter "\\"
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


std::vector<Point> Ifc_clean(std::vector<Point> vertices, std::vector<Face> faces) {
    std::vector<Point> unique_vertices;
    for(int i=0;i<vertices.size();++i)
    {
        if(unique_vertices.empty())
            unique_vertices.emplace_back(vertices[i]);
        else {
            for(int k=0;k<unique_vertices.size();++k)
            {
                if(unique_vertices[k] == vertices[i]) {
                    for(auto & face:faces){
                        for(int & v : face.vertices) {
                            if(v == i)
                            {
                                v = k;
                            }
                        }
                    }
                }
            }
        }
    }
    return unique_vertices;
}


void read_obj(std::string file_in,std::vector<Point> vertices,std::vector<Face> faces)
{
    std::ifstream ifs;

    ifs.open(file_in);
    if(ifs.is_open()==0) {
        std::cout << "can't read input file" << std::endl;
    }
    else if(ifs.is_open()==1) {
        if (ifs.is_open()) {
            std::string line;

            while (getline(ifs, line)) {
                std::istringstream iss(line);
                std::string word;
                iss >> word;
                if(word == "v")
                {
                    std::vector<float> current_vertex;
                    while (iss >> word) current_vertex.push_back(std::stof(word));
                    vertices.emplace_back(Point(current_vertex[0], current_vertex[1], current_vertex[2]));
                }
                if(word == "f")
                {
                    Face face;
                    while(iss>>word){
                        std::string str= split(word, "//")[0];
                        int num=stoi(str);
                        face.vertices.emplace_back(num);
                    }
                    faces.emplace_back(face);
                }
            }
        }
    }
}

int main() {

    std::string data_path="../../data";
    std::string output_path="../../output";

    std::string file_name="/IfcOpenHouse_IFC4.obj";
    std::string file_in = data_path+file_name;

    std::vector<Point> vertices;
    std::vector<Point> uni_vertices;
    std::vector<Face> faces;

    read_obj(file_in,vertices,faces);
    uni_vertices= Ifc_clean(vertices,faces);

    Polyhedron polyhedron;
    Polyhedron_builder<Polyhedron::HalfedgeDS> polyhedron_builder;
    for (auto const &face: ...) {
        polyhedron_builder.faces.emplace_back();
        for (auto const &vertex: face.vertices) {
            polyhedron_builder.vertices.push_back(...);
            polyhedron_builder.faces.back().push_back(...);
        }
    } polyhedron.delegate(polyhedron_builder);
    
    return 0;
}

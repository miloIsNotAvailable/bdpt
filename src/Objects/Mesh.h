#ifndef MESH_H
#define MESH_H

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/vec3.hpp>
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <iostream>
#include <glm/gtc/type_ptr.hpp>
#define TINYOBJLOADER_IMPLEMENTATION
#include "tiny_obj_loader.h"
#include <iostream>
#include <string>

enum class MaterialType {
    Light,
    Diffuse,
    Glossy,
    Mirror,
    Dieletric
};

struct Vertex {
    glm::vec3 position;
    glm::vec3 normal;
    glm::vec3 color;
    MaterialType matType;
    bool light=false;
};

typedef struct {
    std::vector<glm::vec3> vertices;
    std::vector<glm::vec3> colors;
    std::vector<glm::vec3> normals;
    std::vector<float> texcoords;
    std::vector<int> v_indices;
    std::vector<int> vn_indices;
    std::vector<int> vt_indices;
  
    std::vector<tinyobj::material_t> materials;
  
  } MyMesh;  

struct Light {
    public:
    std::vector<glm::vec3> vertices;
    float area;
    tinyobj::shape_t s;

    Light( tinyobj::shape_t &shape, tinyobj::attrib_t& attrib ) {
        computeArea( shape, attrib );
        computeVertices( shape, attrib );
        s = shape;
    }

    private:
    void computeArea(tinyobj::shape_t  &shape, tinyobj::attrib_t& attrib ) {
        
        size_t index_offset = 0.f;
        
        for (size_t f = 0; f < shape.mesh.num_face_vertices.size(); f++) {
            int fv = shape.mesh.num_face_vertices[f];
            if (fv != 3) {
                index_offset += fv; // skip non-triangles
                continue;
            }

            glm::vec3 tri[3];
            for (int v = 0; v < 3; v++) {
                tinyobj::index_t idx = shape.mesh.indices[index_offset + v];
                const float* vp = &attrib.vertices[3 * idx.vertex_index];
                tri[v] = glm::vec3(vp[0], vp[1], vp[2]);
                vertices.push_back(tri[v]);
            }

            glm::vec3 e1 = tri[1] - tri[0];
            glm::vec3 e2 = tri[2] - tri[0];
            area += 0.5f * glm::length(glm::cross(e1, e2));

            index_offset += fv;
        }
    }
    void computeVertices(const tinyobj::shape_t& shape, const tinyobj::attrib_t& attrib) {
        vertices.clear(); // important!
    
        size_t index_offset = 0;
        for (size_t f = 0; f < shape.mesh.num_face_vertices.size(); f++) {
            int fv = shape.mesh.num_face_vertices[f];
            if (fv != 3) {
                index_offset += fv;
                continue; // skip non-triangles
            }
    
            // For each vertex in the face
            for (int v = 0; v < fv; v++) {
                tinyobj::index_t idx = shape.mesh.indices[index_offset + v];
                int v_idx = idx.vertex_index * 3;
    
                glm::vec3 vertex(
                    attrib.vertices[v_idx + 0],
                    attrib.vertices[v_idx + 1],
                    attrib.vertices[v_idx + 2]
                );
    
                vertex /= 200.0f;
                vertex -= glm::vec3(1.0f);
    
                vertices.push_back(vertex);
            }
    
            index_offset += fv;
        }
    }
};

  class Mesh {
    public:
        std::string filename;
        MyMesh mesh;
        std::vector<tinyobj::shape_t> shapes;
        std::vector<tinyobj::material_t> materials;
        std::vector<Light> lights;
        std::vector<Vertex> vertex;

        Mesh(std::string filename) : filename(filename) {
            tinyobj::attrib_t attrib;
            std::string warn, err;
    
            bool ret = tinyobj::LoadObj(&attrib, &shapes, &materials, &warn, &err, filename.c_str(), "../objs/");
    
            if (!ret) {
                std::cerr << "Failed to load OBJ: " << err << std::endl;
                return;
            }
    
            if (!warn.empty()) {
                std::cout << "WARN: " << warn << std::endl;
            }
    
            // Copy materials exactly
            mesh.materials = materials;
    
            // Copy vertices & assign default color (white)
            mesh.vertices.clear();
            mesh.colors.clear();
            for (size_t i = 0; i < attrib.vertices.size(); i += 3) {
                mesh.vertices.emplace_back(attrib.vertices[i], attrib.vertices[i + 1], attrib.vertices[i + 2]);
                mesh.colors.emplace_back(1.0f, 1.0f, 1.0f);  // white color
            }
    
            // Copy normals
            mesh.normals.clear();
            for (size_t i = 0; i < attrib.normals.size(); i += 3) {
                mesh.normals.emplace_back(attrib.normals[i], attrib.normals[i + 1], attrib.normals[i + 2]);
            }
    
            // Copy texcoords (3 floats, z=0)
            mesh.texcoords.clear();
            for (size_t i = 0; i < attrib.texcoords.size(); i += 2) {
                mesh.texcoords.push_back(attrib.texcoords[i]);
                mesh.texcoords.push_back(attrib.texcoords[i + 1]);
                mesh.texcoords.push_back(0.0f);
            }
    
            // Clear indices
            mesh.v_indices.clear();
            mesh.vn_indices.clear();
            mesh.vt_indices.clear();
    
            // FixIndex function is same as your original
            auto FixIndex = [](int idx, int count) {
                if (idx > 0) return idx - 1;
                if (idx < 0) return count + idx;
                return -1;
            };
    
            // Copy indices from shapes
            for (const auto& shape : shapes) {
              for (const auto& index : shape.mesh.indices) {
                  mesh.v_indices.push_back(index.vertex_index);
                  if (index.normal_index >= 0) mesh.vn_indices.push_back(index.normal_index);
                  if (index.texcoord_index >= 0) mesh.vt_indices.push_back(index.texcoord_index);
              }
          }

        for (auto& v : mesh.vertices) {
            v /= 200.0f;
            v -= glm::vec3(1.0f);
        }

          mesh.colors.resize(mesh.vertices.size(), glm::vec3(1.0f)); // default white

          for (size_t s = 0; s < shapes.size(); s++) {
              const auto& shape = shapes[s];
              const auto& mat_ids = shape.mesh.material_ids;
              
              bool isLight = shape.name == "light";
              
              size_t index_offset = 0;
              for (size_t f = 0; f < shape.mesh.num_face_vertices.size(); f++) {
                  int mat_id = mat_ids[f];
                  glm::vec3 color(1.0f);
                  MaterialType mtype = MaterialType::Diffuse;
                  
                  if (mat_id >= 0 && mat_id < materials.size()) {
                    auto& mat = materials[mat_id];
                    bool isDieletric = mat.name == "glass";

                    // std::cout << isDieletric << std::endl; 

                    color = glm::vec3(mat.diffuse[0], mat.diffuse[1], mat.diffuse[2]);

                    if (mat.emission[0] > 0.0f || mat.emission[1] > 0.0f || mat.emission[2] > 0.0f) {
                        mtype = MaterialType::Light;
                    }
                    else if ((mat.specular[0] + mat.specular[1] + mat.specular[2]) > 0.5f) {
                        mtype = MaterialType::Mirror;
                    } 
                    else if( isDieletric ) {
                        mtype = MaterialType::Dieletric;
                    }
                    else {
                        mtype = MaterialType::Diffuse;
                    }
                  }
          
                  for (size_t v = 0; v < 3; v++) {
                    tinyobj::index_t idx = shape.mesh.indices[index_offset + v];
                    
                    int ind = shape.mesh.indices[3 * f + v].vertex_index; 
                    mesh.colors[ind] = color;

                    glm::vec3 pos(
                        attrib.vertices[3 * idx.vertex_index + 0],
                        attrib.vertices[3 * idx.vertex_index + 1],
                        attrib.vertices[3 * idx.vertex_index + 2]
                    );

                    pos /= 200.0f;
                    pos -= glm::vec3(1.0f);
        
                    glm::vec3 normal(0.0f);
                    if (idx.normal_index >= 0) {
                        normal = glm::vec3(
                            attrib.normals[3 * idx.normal_index + 0],
                            attrib.normals[3 * idx.normal_index + 1],
                            attrib.normals[3 * idx.normal_index + 2]
                        );
                    }
        
                    Vertex vert;
                    vert.position = pos;
                    vert.normal   = normal;
                    vert.color    = color;
                    vert.matType  = mtype;
                    vert.light    = isLight;
        
                    // vertex[idx.vertex_index] = vert;
                    vertex.push_back( vert );
                }
        
                index_offset += shape.mesh.num_face_vertices[f];
              }
          }
    
        //   for( tinyobj::shape_t &s : shapes ) {
        //     if( s.name != "light" ) continue;

        //     for (size_t i = 0; i < s.mesh.indices.size(); i++) {
        //       auto idx = s.mesh.indices[i];
  
        //       // Each vertex has 3 floats (x,y,z)
        //       int v_idx = idx.vertex_index * 3;
  
        //       glm::vec3 vertex(
        //           attrib.vertices[v_idx + 0],
        //           attrib.vertices[v_idx + 1],
        //           attrib.vertices[v_idx + 2]
        //       );
  
        //       vertex /= 200.0f;
        //       vertex -= glm::vec3(10.0f);

        //     //   printf( "vec3( %d, %d, %d )\n", vertex.x, vertex.y, vertex.z );
        //     //   lights.push_back(vertex);
        //     }
        //   }

          printf("# of vertices       = %zu\n", mesh.vertices.size());
          printf("# of normals        = %zu\n", mesh.normals.size());
          printf("# of texcoords      = %zu\n", mesh.texcoords.size() / 3);
          printf("# of vertex indices = %zu\n", mesh.v_indices.size());
          printf("# of normal indices = %zu\n", mesh.vn_indices.size());
          printf("# of texcoord indices = %zu\n", mesh.vt_indices.size());


          for( tinyobj::shape_t &shape : shapes ) {
                if( shape.name != "light" ) continue;

                lights.push_back( Light( shape, attrib ) );
                
          } 
        }
    };

#endif
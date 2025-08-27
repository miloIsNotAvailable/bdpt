#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <iostream>
// #include "loadObj.hpp"
#include "intersection.hpp"
// #include "structs.hpp"
#include <optional>
#include "src/Scene/Scene.h"
#include "src/Objects/Model.h"
#include "src/Objects/Shape.h"
#include "src/Camera/Camera.h"
#include <glm/glm.hpp>
#include <cmath>
#include <fstream>
#include <random>
#include <chrono>
#include <thread>

static std::mt19937 rng;

// initialize once at program start
void initRandom() {
    // seed with high-resolution clock
    uint64_t seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    rng.seed(static_cast<unsigned>(seed));
}

float randFloat(float a, float b) {
    std::uniform_real_distribution<float> dist(a, b);
    return dist(rng);
}

int randInt(int a, int b) {
    std::uniform_int_distribution<int> dist(a, b);
    return dist(rng);
}


glm::vec3 randomPointOnHemisphere(const glm::vec3& normal, float r1, float r2) {
    // Sample in local space (hemisphere around +Z)
    float phi = 2.0f * glm::pi<float>() * r1;
    float cosTheta = std::sqrt(1.0f - r2);
    float sinTheta = std::sqrt(r2);

    glm::vec3 localDir(
        std::cos(phi) * sinTheta,
        std::sin(phi) * sinTheta,
        cosTheta
    );

    // Build an orthonormal basis around the normal
    glm::vec3 up = std::abs(normal.z) < 0.999f ? glm::vec3(0.0f, 0.0f, 1.0f)
                                               : glm::vec3(1.0f, 0.0f, 0.0f);

    glm::vec3 tangent = glm::normalize(glm::cross(up, normal));
    glm::vec3 bitangent = glm::cross(normal, tangent);

    // Transform local direction to world space
    return glm::normalize(
        tangent * localDir.x +
        bitangent * localDir.y +
        normal * localDir.z
    );
}

struct Pixel {
    glm::vec2 p;
    glm::vec3 L;

    Pixel( glm::vec2 p, glm::vec3 L ): p(p), L(L) {}
    Pixel( ): p(glm::vec2(0.f)), L(glm::vec3( 0.f )) {}
};

void savePPM(const std::string &filename,
    const std::vector<Pixel> &framebuffer,
    int width, int height)
{
    std::ofstream ofs(filename, std::ios::binary);
    ofs << "P6\n" << width << " " << height << "\n255\n";

    for (int i = 0; i < width * height; ++i) {
        glm::vec3 c = glm::clamp(framebuffer[i].L, 0.0f, 1.0f);
        unsigned char r = static_cast<unsigned char>(c.r * 255.0f);
        unsigned char g = static_cast<unsigned char>(c.g * 255.0f);
        unsigned char b = static_cast<unsigned char>(c.b * 255.0f);
        ofs.write(reinterpret_cast<char*>(&r), 1);
        ofs.write(reinterpret_cast<char*>(&g), 1);
        ofs.write(reinterpret_cast<char*>(&b), 1);
    }
}

struct PathVertex {
    glm::vec3 vec;
    glm::vec3 dir;
    glm::vec3 nor;
    glm::vec3 beta;
    glm::vec3 alpha = glm::vec3( 1.f );
    float pdf=1.f;
    float p=1.f;
    float pA=1.f;
    glm::vec3 color;
    float pdfFwd=1.; 
    float pdfRev=1.; 
    PathVertex* prev;
    bool isLight=false;
    MaterialType mtype = MaterialType::Diffuse;

    PathVertex(glm::vec3 v, glm::vec3 dir, glm::vec3 nor, glm::vec3 beta = glm::vec3(1.f))
    : vec(v), dir(dir), nor(nor), beta(beta) {}

    PathVertex(Hit &hitPoint, const glm::vec3 &betaIn = glm::vec3(1.f), float eps = 1e-4f) {
        vec = hitPoint.hitPoint + hitPoint.hitNormal * eps;
        nor = hitPoint.hitNormal;
        dir = randomPointOnHemisphere(hitPoint.hitNormal, glm::linearRand(0.f, 1.f), glm::linearRand(0.f, 1.f));
        beta = betaIn;
    }
};

void BSDF( MaterialType &m, glm::vec3 &wi, 
           glm::vec3 &n, glm::vec3 &newDir, 
           float &pdf, float &cosTheta, glm::vec3 &throughput, glm::vec3 &color ) {
    switch( m ) {
        // case MaterialType::Diffuse: {}

        case MaterialType::Dieletric: {
            // float ior = 1.5;
            // float F0 = pow((1.0 - ior) / (1.0 + ior), 2.0);

            float n1 = 1.;
            float n2 = 1.5;

            float cos_theta =glm::dot(-wi, n);

            if (cos_theta < 0.) {
                // Swap the indices of refraction for exiting the medium
                // vec3 nor_ = -nor;
                n = -n;
                cos_theta =glm::dot(-wi, n);
                // float temp = n1;
                n1 = n2;
                // n2 = temp;
                // break;

            }

            float eta = n1 / n2;
            float k = 1.0 - eta * eta * (1.0 - cos_theta * cos_theta);

            float reflect_prob, R0;
            glm::vec3 reflect_dir, refract_dir;

            if (k < 0.0) {
                reflect_prob = 1.0;
                // refract_dir = vec3(.0);
                reflect_dir = glm::normalize(glm::reflect(wi, n));

            } else {
                refract_dir = glm::refract( wi, n, eta );
                R0 = ((n1 - n2) / (n1 + n2)) * ((n1 - n2) / (n1 + n2));
                reflect_prob = R0 + (1.0 - R0) * pow(1.0 - cos_theta, 5.0);
                reflect_dir = glm::normalize(glm::reflect(wi, n));
            }

            if (randFloat(0., 1.) < reflect_prob) {
                newDir = reflect_dir ;
            } else {
                newDir = refract_dir;
            }
            cosTheta = cos_theta;
            pdf = 1.;

            throughput *= color;

            break;
        }

        default: {
            glm::vec3 dir = randomPointOnHemisphere(n,
                randFloat(0.f,1.f),
                randFloat(0.f,1.f));

            cosTheta = std::max(0.f, glm::dot(n, dir));
            pdf = cosTheta / glm::pi<float>();
            newDir = dir;
            // if (pdf < 1e-8f || cosTheta <= 0.f) break;
            glm::vec3 albedo(color);
            glm::vec3 f = color / glm::pi<float>();              // Lambertian f_s
    
            throughput *= f * (cosTheta / pdf);
            break;
        }
    };
}

std::vector<std::vector<glm::vec3>> cameraPath;
std::vector<std::vector<glm::vec3>> lightPath;

PathVertex generateCameraPoint(Camera& c, const glm::vec2& pixel_i) {
    glm::vec3 v = c.screenToWorldSpace( pixel_i );
    glm::vec3 p = c.c_sensor + v;
    glm::vec3 rayDir = glm::normalize( p - c.center);
    
    glm::vec3 F = c.getFocalPlane( c.focalDist );
    float t = c.getFocalPointFactor( F, rayDir );

    glm::vec3 focalPoint = c.center + t * rayDir;

    glm::vec2 rndAperturePoint = c.sampleRndPointOnDisk( c.apertureSize );
    glm::vec3 rndAperturePointWorld = c.screenToWorldSpace( rndAperturePoint );

    glm::vec3 pointOnAperture = c.center + rndAperturePointWorld;

    glm::vec3 pRayDir = glm::normalize( focalPoint - pointOnAperture );
    
    float r2   = glm::length2(focalPoint - pointOnAperture);
    float cosL = std::abs(glm::dot(c.eye, pRayDir));    // lens plane normal

    float pA_lens = 1.f / (glm::pi<float>() );
    float p_omega = pA_lens * (r2 / std::max(1e-8f, cosL));  // area → solid-angle
    float pdf_pixel = 1./ (400. * 400.); 

    PathVertex z0(pointOnAperture, pRayDir, c.eye);
    z0.pdfFwd = p_omega; 
    z0.pdfRev = pdf_pixel;       
    z0.p      = p_omega;   
    float pdf_ray = pdf_pixel * p_omega;
    z0.alpha  = glm::vec3(1.f);
    return z0;
}

PathVertex generateRandomLightPoint( Model &m, float &pdf ) {

    int rndlightInd = randFloat( 0.f, float(m.mesh.lights.size() - 1) );
    Light l = m.mesh.lights[ rndlightInd ];

    int triCount = l.vertices.size() / 3;
    int i = randInt(0, triCount - 1);

    glm::vec3 v0 = l.vertices[(3 * i + 0)];
    glm::vec3 v1 = l.vertices[(3 * i + 1)];
    glm::vec3 v2 = l.vertices[(3 * i + 2)];

    // std::cout << i << std::endl;

    glm::vec3 lightNormal = glm::normalize(glm::cross( v0 - v2, v0 - v1 ));
    float triArea = 0.5f * glm::length(glm::cross(v1 - v0, v2 - v0));

    float pdfA = 1.0f / l.area;
    float lightPickProb = 1./m.mesh.lights.size();

    float r1 = randFloat( 0.f, 1.f );
    float r2 = randFloat( 0.f, 1.f );

    float u = 1.f - sqrt(r1);
    float v = r2 * sqrt(r1);
    
    // float alpha = 1 - u - v;        

    glm::vec3 p = (1 - u - v) * v0 + u * v1 + v * v2;

    glm::vec3 down = glm::vec3( 0., -1., 0. );
    bool isDown = glm::dot( down, lightNormal ) < 0.;

    // printf( "normal: %f, %f, %f\n", lightNormal.x, lightNormal.y, lightNormal.z );

    if( isDown ) lightNormal *= -1.;

    glm::vec3 dir = randomPointOnHemisphere( lightNormal, 
                                             glm::linearRand(0.f, 1.f), 
                                             glm::linearRand(0.f, 1.f) );
    

    float cosTheta = std::max( 0.f, glm::dot( lightNormal, dir ) );
    float pdfW = cosTheta / glm::pi<float>();

    glm::vec3 Le = glm::vec3( 20.f );
    glm::vec3 alpha = Le / (pdfA * pdfW * lightPickProb);
    // glm::vec3 alpha = Le;

    PathVertex pv( p, dir, lightNormal );

    pv.alpha   = Le;
    pv.pdfFwd  = pdfW;
    pv.pdfRev  = pdfA * lightPickProb;
    pv.p       = pdfW * pdfA * lightPickProb;
    pv.color   = glm::vec3( 20.f );

    // pdf = pdfW * pdfA * lightPickProb;
    
    return pv;
}

void generateCameraSubPath(Model &m, Camera &c, glm::vec2 &pixel_i,
    std::vector<PathVertex> &path, int depth) 
{
    PathVertex z0 = generateCameraPoint(c, pixel_i);
    
    path.push_back( z0 );
    
    glm::vec3 throughput = z0.alpha;
    PathVertex s = z0;
    PathVertex prev = s;

    for( int i = 0; i < depth; i ++ ) {
        auto hitOpt = m.isect(s.vec, s.dir);
        if (!hitOpt.has_value()) break;
        // if (hitOpt.value().mtype == MaterialType::Light) break; 

        Hit hit = hitOpt.value();
        glm::vec3 pos = hit.hitPoint;
        glm::vec3 n   = hit.hitNormal;

        // bool isD = hit.mtype == MaterialType::Dieletric;
        // std::cout << isD;

        float pdf, cosTheta;
        glm::vec3 newDir;
        BSDF( hit.mtype, s.dir, n, newDir, pdf, cosTheta, throughput, hit.color );

        // glm::vec3 newDir = randomPointOnHemisphere(n,
        //     randFloat(0.f,1.f),
        //     randFloat(0.f,1.f));

        // float cosTheta = std::max(0.f, glm::dot(n, newDir));
        // float pdf = cosTheta / glm::pi<float>();
        if (pdf < 1e-8f || cosTheta <= 0.f) break;
        
        float pdfFwd=1.f;
        float pdfRev=1.f;
        if( i < 1 ) {
            pdfFwd = z0.pdfFwd;
            pdfRev = z0.pdfRev;
        } else {
            glm::vec3 dir = s.vec - prev.vec;
            glm::vec3 w   = glm::normalize(dir);
            
            float dist2 = glm::dot(dir, dir);
            float cosX = std::max(0.f, glm::dot(prev.nor, w));
            float cosY = std::max(0.f, glm::dot(s.nor, -w));
            float G = (cosX * cosY) / dist2;

            float pdf = cosX / glm::pi<float>();
    
            pdfRev = G * pdf;

            glm::vec3 dirNew = pos - s.vec;
            glm::vec3 wNew = glm::normalize(dirNew);
            float dist2New = glm::dot(dirNew, dirNew);

            float cosXNew = std::max(0.f, glm::dot(s.nor, wNew));
            float cosYNew = std::max(0.f, glm::dot(n, -wNew));
            float GNew = (cosXNew * cosYNew) / dist2New;

            float pdfNew = cosXNew / glm::pi<float>();
    
            pdfFwd = GNew * pdfNew;
             
        }

        // printf( "pdf rev: %f, pdf fwd: %f\n", pdfRev, pdfFwd );

        // float cosThetaOld = std::max(0.f, glm::dot(n, s.dir));
        // float pdfOld = cosThetaOld / glm::pi<float>();              // solid-angle
        // if (pdfOld < 1e-8f || cosThetaOld <= 0.f) break;

        glm::vec3 albedo(hit.color);
        // glm::vec3 f = hit.color / glm::pi<float>();              // Lambertian f_s

        // throughput *= f * (cosTheta / pdf);
        // throughput *= f;

        PathVertex pv(pos, newDir, n);
        pv.alpha  = throughput;
        pv.color  = albedo;
        pv.pdfFwd = pdfFwd;
        pv.pdfRev = pdfRev;
        pv.p      = pdf;
        pv.mtype  = hit.mtype;
        pv.isLight= hit.mtype == MaterialType::Light;

        // float pCont = std::min(0.95f, glm::compMax(throughput));
        // if (randFloat(0.f,1.f) > pCont) break;
        // throughput /= pCont;

        if( pv.isLight ) {
            pv.color = glm::vec3( 20.f );
            path.push_back(pv);
            break;
        }


        path.push_back(pv);
        // if (hit.mtype == MaterialType::Light) break;

        prev = s;
        s = pv;
    }

    // // We / P_A(z0), assume We as 1. for now
    // s.alpha = glm::vec3(1.f);
    // // P_A(z0) on thin lens
    // s.p = 1.f/(glm::pi<float>() * c.apertureSize * c.apertureSize); 
    // path.push_back(s);
    
    // glm::vec3 f( 1.f );
    // PathVertex prev = s;

    // for (int i = 1; i < depth; ++i) {
    //     auto hitPoint = m.isectRay(s.vec, s.dir);
    //     if (!hitPoint.has_value()) break;

    //     Hit &hit = hitPoint.value();
    //     glm::vec3 pos = hit.hitPoint + hit.hitNormal * 1e-4f;
    //     glm::vec3 normal = hit.hitNormal;

    //     glm::vec3 newDir = randomPointOnHemisphere(normal,
    //                 glm::linearRand(0.f, 1.f),
    //                 glm::linearRand(0.f, 1.f));

    //     float cosTheta = std::max(0.f, glm::dot(normal, newDir));
    //     // if (cosTheta <= 0.f) break;

    //     float pdf = cosTheta / glm::pi<float>();
        
    //     // last alpha * f_s * cosTheta / pA
    //     glm::vec3 newAlpha = prev.alpha * f * cosTheta / prev.p;

    //     glm::vec3 albedo( hit.color );
    //     f = albedo / glm::pi<float>();
        
    //     PathVertex pv(pos, newDir, normal);
    //     pv.p = pdf;
    //     pv.pA = pdf;
    //     pv.color=albedo;
    //     pv.alpha = newAlpha;

    //     path.push_back(pv);

    //     if( hitPoint.value().mtype == MaterialType::Light ) break;

    //     // prev = s;
    //     s = pv; 
    // }
}

void generateLightSubPath(Model &m, std::vector<PathVertex> &path, int depth) 
{
    float pdf;
    PathVertex y0 = generateRandomLightPoint(m, pdf);
    path.push_back( y0 );
    
    glm::vec3 throughput = y0.alpha;
    // printf( "%f, %f, %f\n", throughput.x, throughput.y, throughput.z );
    PathVertex s = y0;
    PathVertex prev = s;

    for( int i = 0; i < depth; i ++ ) {
        auto hitOpt = m.isect(s.vec, s.dir);
        if (!hitOpt.has_value()) break;
        // if (hitOpt.value().mtype == MaterialType::Light) break; 

        Hit hit = hitOpt.value();
        glm::vec3 pos = hit.hitPoint;
        glm::vec3 n   = hit.hitNormal;

        // if (hit.mtype == MaterialType::Light) break;

        float pdf, cosTheta;
        glm::vec3 newDir;
        BSDF( hit.mtype, s.dir, n, newDir, pdf, cosTheta, throughput, hit.color );

        // glm::vec3 newDir = randomPointOnHemisphere(n,
        //     randFloat(0.f,1.f),
        //     randFloat(0.f,1.f));

        // float cosTheta = std::max(0.f, glm::dot(n, newDir));
        // float pdf = cosTheta / glm::pi<float>();
        // if (pdf < 1e-8f || cosTheta <= 0.f) break;
        
        float pdfFwd = 1.f;
        float pdfRev = 1.f;
        if( i < 1 ) {
            pdfFwd = y0.pdfFwd;
            pdfRev = y0.pdfRev;
        } else {
            glm::vec3 dir = s.vec - prev.vec;
            glm::vec3 w   = glm::normalize(dir);
            
            float dist2     = glm::dot(dir, dir);
            float cosX = std::max(0.f, glm::dot(prev.nor, w));
            float cosY = std::max(0.f, glm::dot(s.nor, -w));
            float G = (cosX * cosY) / dist2;

            float pdf = cosX / glm::pi<float>();
    
            pdfRev = G * pdf;

            glm::vec3 dirNew = pos - s.vec;
            glm::vec3 wNew = glm::normalize(dirNew);
            float dist2New = glm::dot(dirNew, dirNew);

            float cosXNew = std::max(0.f, glm::dot(s.nor, wNew));
            float cosYNew = std::max(0.f, glm::dot(n, -wNew));
            float GNew = (cosXNew * cosYNew) / dist2New;

            float pdfNew = cosXNew / glm::pi<float>();
    
            pdfFwd = GNew * pdfNew;
        }


        // float cosThetaOld = std::max(0.f, glm::dot(n, s.dir));
        // float pdfOld = cosThetaOld / glm::pi<float>();              // solid-angle
        // if (pdfOld < 1e-8f || cosThetaOld <= 0.f) break;

        glm::vec3 albedo(hit.color);
        // glm::vec3 f = albedo / glm::pi<float>();              // Lambertian f_s

        // throughput *= f * (cosTheta / pdf);
        // throughput *= f;

        PathVertex pv(pos, newDir, n);
        pv.alpha  = throughput;
        pv.color  = albedo;
        pv.pdfFwd = pdfFwd;
        pv.pdfRev = pdfRev;
        pv.p      = pdf;
        pv.mtype  = hit.mtype;
        pv.isLight= hit.mtype == MaterialType::Light;

        if( pv.isLight ) {
            pv.color = glm::vec3( 20.f );
            path.push_back(pv);
            break;
        }

        // printf( "%f, %f, %f\n", hit.color.x, hit.color.y, hit.color.z );
        // printf( "%f, %f, %f\n\n", throughput.x, throughput.y, throughput.z );

        // float pCont = std::min(0.95f, glm::compMax(throughput));
        // if (randFloat(0.f,1.f) > pCont) break;
        // throughput /= pCont;

        path.push_back(pv);

        prev = s;
        s = pv;
    }


    // float pdf;
    // PathVertex s = generateRandomLightPoint(m, pdf);
    
    // // We / P_A(z0), assume We as 1. for now
    // s.alpha = glm::vec3(1.f);
    // // P_A(z0) on thin lens
    // s.p = 1.f/(glm::pi<float>()); 
    // path.push_back(s);

    // glm::vec3 f( 1.f );
    // PathVertex prev = s;

    // for (int i = 1; i < depth; ++i) {
    //     auto hitPoint = m.isectRay(s.vec, s.dir);
    //     if (!hitPoint.has_value()) break;
    //     // if( hitPoint.value().mtype == MaterialType::Light ) break;

    //     Hit &hit = hitPoint.value();
    //     glm::vec3 pos = hit.hitPoint + hit.hitNormal * 1e-4f;
    //     glm::vec3 normal = hit.hitNormal;

    //     // Diffuse sampling
    //     glm::vec3 newDir = randomPointOnHemisphere(normal,
    //                 glm::linearRand(0.f, 1.f),
    //                 glm::linearRand(0.f, 1.f));

    //                 float cosTheta = std::max(0.f, glm::dot(normal, newDir));
    //     // if (cosTheta <= 0.f) break;

    //     // last alpha * f_s * cosTheta / pA
    //     glm::vec3 newAlpha = prev.alpha * f * cosTheta / prev.p;

    //     glm::vec3 albedo( hit.color );
    //     float PI = 3.1415926f;
    //     f = albedo / PI;
    //     float pdf = cosTheta / PI;

    //     PathVertex pv(pos, newDir, normal);
    //     pv.p = pdf;
    //     pv.pA = pdf;
    //     pv.color=albedo;
    //     pv.alpha = newAlpha;

    //     path.push_back(pv);

    //     prev = s;
    //     s = pv; 
    // }
}

float MISWeight(std::vector<PathVertex> &l, std::vector<PathVertex> &e, int s, int t) {
    std::vector<float> pdfFwd(s+t);
    std::vector<float> pdfRev(s+t);

    // copy forward/backward pdfs into a single path
    for(int i=0; i<s; i++) {
        pdfFwd[i] = l[i].pdfFwd;
        pdfRev[i] = l[i].pdfRev;

        // printf( "pdf rev: %f, pdf fwd: %f\n", l[i].pdfRev, pdfFwd[i] );
    }
    for(int i=0; i<t; i++) {
        pdfFwd[s+i] = e[i].pdfFwd;
        pdfRev[s+i] = e[i].pdfRev;
    }

    // compute the balance heuristic
    float sum = 1.f;  // start with chosen strategy
    float prod = 1.f;

    // iterate over all alternative strategies
    for(int i=0; i<s+t-1; i++) {
        prod *= pdfRev[i] / pdfFwd[i];
        sum += prod;
        // printf( "pdf rev: %f, pdf fwd: %f\n", pdfRev[i], pdfFwd[i] );
    }

    // std::cout << sum << "\n";
    if( sum > 0.f ) 
        return 1.f / sum;
    else 
        return 1.;
}

void connectPath( Model &m, Camera &c, glm::vec2 &p, glm::vec3 &L ) {
    std::vector<PathVertex> ePath;
    std::vector<PathVertex> lPath;

    generateCameraSubPath( m, c, p, ePath, 5 );
    generateLightSubPath( m, lPath, 5 );

    for (int t = 0; t < ePath.size(); t++) {
        PathVertex& x = ePath[t];
        for (int s = 0; s < lPath.size(); s++) {
            PathVertex& y = lPath[s];
            
            // printf( "%f, %f\n", x.pdfFwd, x.pdfRev );

            glm::vec3 dir = (y.vec) - (x.vec);
            float dist2   = glm::dot(dir, dir);
            dir = glm::normalize(dir);
    
            // printf( "%f\n", dist2 );

            float tmin = 0.f;
            float tmax = std::sqrt( dist2 );

            // visibility check
            if (m.isect(x.vec, dir, &tmin, &tmax)) continue;
    
            glm::vec3 c;
            if( s == 0 && t > 0 ) {
                if( y.isLight ) {
                    c = glm::vec3( y.color );
                }
                //     c = glm::vec3( 0.f );
            } else if( t > 0 && s > 0 ) {

                // evaluate BSDFs
                glm::vec3 f_x = x.color / glm::pi<float>();
                glm::vec3 f_y = y.color / glm::pi<float>();

                // geometry term
                float cosX = std::max(0.f, glm::dot(x.nor, dir));
                float cosY = std::max(0.f, glm::dot(y.nor, -dir));
                float G = (cosX * cosY) / dist2;
        
                // contribution
                c = f_x * G * f_y;

                glm::vec3 contrib = x.alpha * c * y.alpha;
        
                // MIS weight (placeholder = 1)
                float w = MISWeight( lPath, ePath, s, t ); // TODO: compute with balance heuristic
                // printf( "%f\n", w );
                L += w * contrib;
                // printf( "w: %f, L: %f, %f, %f\n", w, y.alpha.x, y.alpha.y, y.alpha.z );
                
                // glm::vec3 f_x = evalBSDF(x, dir);        // f(x→y)
                // glm::vec3 f_y = evalBSDF(y, -dir);       // f(y→x), or emission if y is light
        
            }
            // printf( "%f, %f, %f\n\n", y.vec.x, y.vec.y, y.vec.z );

        }
    }
}

int main() {

    initRandom();

    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    int w = 800;
    int h = 800;
    int DEPTH = 5;

    GLFWwindow* win = glfwCreateWindow(w, h, "Quad", nullptr, nullptr);
    if (!win) return -1;
    glfwMakeContextCurrent(win);

    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
        std::cerr << "Failed to load GLAD\n";
        return -1;
    }

    
    Mesh mesh = Mesh( "../objs/cornell_box.obj" );
    
    Scene s = Scene( w, h );
    Model m( mesh );

    // for( auto &t : m.bvh->triangles) {
    //     printf( "%f, %f, %f\n", mesh.vertex[ t.idx.x ].position.x, mesh.vertex[ t.idx.x ].position.y, mesh.vertex[ t.idx.x ].position.z );
    //     printf( "%f, %f, %f\n", mesh.vertex[ t.idx.y ].position.x, mesh.vertex[ t.idx.y ].position.y, mesh.vertex[ t.idx.y ].position.z );
    //     printf( "%f, %f, %f\n\n", mesh.vertex[ t.idx.z ].position.x, mesh.vertex[ t.idx.z ].position.y, mesh.vertex[ t.idx.z ].position.z );
    // }

    glm::vec3 center( .5, 1.0, -6. );
    glm::vec3 eye( .3, .3, 0. );
    float apertureSize = .1;
    float focalDist = 5.5;
    float c_sensor_dist = 2.f;
    Camera c( center, eye, apertureSize, focalDist, c_sensor_dist );
    c.setCameraBasis();

    Shape shape( GL_POINTS );
    shape.beginShape();
    for( int i = 0; i < 100; i ++ ) {
        glm::vec2 v = c.sampleRndPointOnDisk( apertureSize );
        glm::vec3 vWorld = c.screenToWorldSpace( v );

        shape.vertex( center + vWorld );
        shape.color( glm::vec3( 1., 1., 0. ) );
    }
    shape.endShape();

    glm::vec3 e1( 1., 1., 0. );
    glm::vec3 e2( 1., -1., 0. );
    glm::vec3 e3( -1., 1., 0. );
    glm::vec3 e4( -1., -1., 0. );

    glm::vec3 p1( c.focalPlane + c.screenToWorldSpace( e1 ) );
    glm::vec3 p2( c.focalPlane + c.screenToWorldSpace( e2 ) );
    glm::vec3 p3( c.focalPlane + c.screenToWorldSpace( e3 ) );
    glm::vec3 p4( c.focalPlane + c.screenToWorldSpace( e4 ) );

    Shape focalPlane( GL_TRIANGLES );
    focalPlane.beginShape();
    focalPlane.vertex( p1 ); focalPlane.color( glm::vec3( 1. ) );
    focalPlane.vertex( p2 ); focalPlane.color( glm::vec3( 1. ) );
    focalPlane.vertex( p3 ); focalPlane.color( glm::vec3( 1. ) );

    focalPlane.vertex( p2 ); focalPlane.color( glm::vec3( 1. ) );
    focalPlane.vertex( p3 ); focalPlane.color( glm::vec3( 1. ) );
    focalPlane.vertex( p4 ); focalPlane.color( glm::vec3( 1. ) );
    focalPlane.endShape();

    std::vector<glm::vec2> pixel;
    int GRID_SIZE = 8;
    int SIZE = GRID_SIZE * GRID_SIZE;
    float width = 500.;
    float height = 500.;

    for( int i = 0; i < SIZE; i ++ ) {
    
        int x = (i % GRID_SIZE) * std::floor(width/GRID_SIZE) - width/2;
        int y = ((int)(std::floor( i / GRID_SIZE )) % GRID_SIZE) * std::floor(height/GRID_SIZE) - height/2;
        // console.log( x, y ,z )
    
        pixel.push_back(glm::vec2( ((float)x / (float)width), ((float)y / (float)height) ));
    }

    Shape pixels(GL_POINTS);
    pixels.beginShape();
    for( int i = 0; i < SIZE; i ++ ) {
        glm::vec3 v = c.screenToWorldSpace( pixel[i] );
        pixels.vertex( c.c_sensor + v );
        pixels.color( glm::vec3( 0., 1., 0. ) );
    }
    pixels.endShape();

    std::vector<glm::vec3> pointsOnLight;
    Shape lRays(GL_LINES);
    Shape rays(GL_LINES);
    Shape connected( GL_LINES );

    for( int i = 0; i < SIZE; i ++ ) {
        std::vector<glm::vec3> cPath;
        glm::vec3 v = c.screenToWorldSpace( pixel[i] );
        glm::vec3 p = c.c_sensor + v;
        glm::vec3 rayDir = glm::normalize( p - c.center);
        
        glm::vec3 F = c.getFocalPlane( c.focalDist );
        float t = c.getFocalPointFactor( F, rayDir );

        glm::vec3 focalPoint = c.center + t * rayDir;

        glm::vec2 rndAperturePoint = c.sampleRndPointOnDisk( apertureSize );
        glm::vec3 rndAperturePointWorld = c.screenToWorldSpace( rndAperturePoint );

        glm::vec3 pointOnAperture = c.center + rndAperturePointWorld;

        glm::vec3 pRayDir = glm::normalize( focalPoint - pointOnAperture );
        
        cPath.push_back( pointOnAperture );

        for( int i = 0; i < DEPTH; i ++ ) {
            
            auto hitPoint = m.isectRay( pointOnAperture, pRayDir );
            
            if( !hitPoint.has_value() ) break;

            // glm::vec3 hitPoint = pointOnAperture + (closestHit.value()) * pRayDir;

            // rays.vertex( pointOnAperture );
            // rays.color( glm::vec3( 0., 1., 1. ) );
            // rays.vertex( hitPoint );
            // rays.color( glm::vec3( 0., 1., 1. ) );

            cPath.push_back( hitPoint.value().hitPoint );

            // pointOnAperture = hitPoint + hitNormal * 0.00001f;
            pointOnAperture = hitPoint.value().hitPoint;
            pRayDir = randomPointOnHemisphere( hitPoint.value().hitNormal, glm::linearRand( 0., 1. ), glm::linearRand( 0., 1. ) );
        }

        cameraPath.push_back( cPath );
    }

    rays.beginShape();
    lRays.beginShape();
    connected.beginShape();

    for( int i = 0; i < 3; i ++ ) {
        glm::vec2 p = pixel[i];
        std::vector<PathVertex> cPath;
        std::vector<PathVertex> lPath;
        generateCameraSubPath( m, c, p, cPath, 5 );
        generateLightSubPath( m, lPath, 5 );
    
        for( int i = 1; i < cPath.size(); i++ ) {
            glm::vec3 prev = cPath[i - 1].vec;
            glm::vec3 curr = cPath[i].vec;
            
            rays.vertex(prev);
            rays.color(glm::vec3(0., 1., 1.));
            rays.vertex(curr);
            rays.color(glm::vec3(0., 1., 1.));
        }

        for( int i = 1; i < lPath.size(); i++ ) {
            glm::vec3 prev = lPath[i - 1].vec;
            glm::vec3 curr = lPath[i].vec;
            
            lRays.vertex(prev);
            lRays.color(glm::vec3(1., 0., 1.));
            lRays.vertex(curr);
            lRays.color(glm::vec3(1., 0., 1.));
        }

        for (int t = 0; t < cPath.size(); t++) {
            PathVertex& x = cPath[t];
            for (int s = 0; s < lPath.size(); s++) {

                PathVertex& y = lPath[s];
                    
                glm::vec3 dir = (y.vec) - (x.vec);
                float dist2   = glm::dot(dir, dir);
                dir = glm::normalize(dir);
        
                float tmin = 1e-4f;
                float tmax = std::sqrt( dist2 ) - tmin;

                if (m.isect(x.vec, dir, &tmin, &tmax)) continue;

                if ( s == 0 && t > 0 ){
                    connected.vertex( x.vec ); connected.color( glm::vec3( 1., 0., 1. ) );
                    connected.vertex( y.vec ); connected.color( glm::vec3( 1., 0., 1. ) );

                } else if( s > 0 && t > 0 ) {

                    connected.vertex( x.vec ); connected.color( glm::vec3( 1., 0., 0. ) );
                    connected.vertex( y.vec ); connected.color( glm::vec3( 1., 0., 0. ) );

                    // float w = MISWeight( lPath, cPath, s, t );
                    // printf( "%f\n", w );
                }
            }
        }
    }
    rays.endShape();
    lRays.endShape();
    connected.endShape();

    while (!glfwWindowShouldClose(win)) {
        glfwPollEvents();
        s.orbitControls( win );

        glClearColor(0.1f, 0.1f, 0.15f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glEnable(GL_DEPTH_TEST);

        s.draw();
        m.draw();
        glPointSize(5.0f);
        shape.draw();
        pixels.draw();
        // focalPlane.draw();
        // rays.draw();
        lRays.draw();
        // lightPoints.draw();
        // connected.draw();

        glfwSwapBuffers(win);
    }
    glfwDestroyWindow(win);
    glfwTerminate();

    std::vector<Pixel> imgPixels;
    GRID_SIZE = 400;
     SIZE = GRID_SIZE * GRID_SIZE;
     width = 400.;
     height = 400.;

    for( int i = 0; i < SIZE; i ++ ) {
    
        int x = (i % GRID_SIZE) * std::floor(width/GRID_SIZE) - width/2;
        int y = ((int)(std::floor( i / GRID_SIZE )) % GRID_SIZE) * std::floor(height/GRID_SIZE) - height/2;
        // console.log( x, y ,z )
    
        Pixel p( glm::vec2( (float)x / (float)width, (float)y / (float)height ), 
                 glm::vec3( 0.f ) );

        imgPixels.push_back( p );
    }


    int numThreads = std::thread::hardware_concurrency();
    std::vector<std::thread> workers;

    int chunkSize = imgPixels.size() / numThreads;

    auto connect = [&](std::vector<Pixel>::iterator startIt, std::vector<Pixel>::iterator endIt) {
        for (auto it = startIt; it != endIt; ++it) {
            Pixel &p = *it;
            int iter = 120;
            for (int i = 0; i < iter; i++) {
                connectPath(m, c, p.p, p.L);
            }
    
            glm::vec3 col = p.L / float(iter);
            col = col / (col + glm::vec3(1.0));
            col = glm::pow(col, glm::vec3(1.0 / 2.2));
            p.L = col;
        }
    };
    
    for (int t = 0; t < numThreads; t++) {
        auto startIt = imgPixels.begin() + t * chunkSize;
        auto endIt = (t == numThreads - 1) ? imgPixels.end() : startIt + chunkSize;
    
        workers.emplace_back(connect, startIt, endIt);
    }
    
    for (auto &worker : workers) {
        worker.join();
    }

    // std::vector<glm::vec3> fb;
    // for( Pixel &p : imgPixels ) {
    //     // glm::vec3 L(0.0f);

    //     int iter = 10;
    //     for( int i = 0; i < iter; i ++ ) {
    //         connectPath( m, c, p.p, p.L );
    //     }

    //     glm::vec3 col = p.L / float(iter);

    //     // float exposure = 2.5f;
    //     // col = 1.0f - glm::exp(-col * exposure);

    //     col = col / (col + glm::vec3(1.0));
    //     col = glm::pow(col, glm::vec3(1.0 / 2.2));

    //     p.L = col;

    //     // fb.push_back(col);    
    // }

    savePPM( "hi.ppm", imgPixels, width, height );

    // glDeleteVertexArrays(1, &VAO);
    // glDeleteBuffers(1, &VBO);
    m.~Model();
    // glDeleteProgram(shader);
    s.~Scene();
    return 0;
}
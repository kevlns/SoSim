//
// Created by ADMIN on 2024/3/13.
//

#include <fstream>

#include "libs/ModelL/model_exporter.hpp"

namespace SoSim {

    void ModelExporter::exportVecSetAsPly(const std::string &dir,
                                          const std::string &file_name,
                                          const std::vector<Vec3f> &pos) {
#ifdef WIN32
        std::ofstream ofs(dir + '\\' + file_name + ".ply");
#else
        std::ofstream ofs(dir + '/' + file_name + ".ply");
#endif

        ofs << "ply\n";
        ofs << "format ascii 1.0\n";
        ofs << "element vertex " << pos.size() << "\n";
        ofs << "property float x\n";
        ofs << "property float y\n";
        ofs << "property float z\n";
        ofs << "end_header\n";

        for (auto &po: pos) {
            ofs << po.x << " " << po.y << " " << po.z << "\n";
        }

        ofs.close();
    }

    void ModelExporter::exportVecSetAsPly(const std::string &dir,
                                          const std::string &file_name,
                                          const std::vector<Vec3f> &pos,
                                          const std::vector<Vec3f> &color) {
#ifdef WIN32
        std::ofstream ofs(dir + '\\' + file_name + ".ply");
#else
        std::ofstream ofs(dir + '/' + file_name + ".ply");
#endif

        ofs << "ply\n";
        ofs << "format ascii 1.0\n";
        ofs << "element vertex " << pos.size() << "\n";
        ofs << "property float x\n";
        ofs << "property float y\n";
        ofs << "property float z\n";
        ofs << "property uchar red" << std::endl;
        ofs << "property uchar green" << std::endl;
        ofs << "property uchar blue" << std::endl;
        ofs << "end_header\n";

        for (int i = 0; i < pos.size(); ++i) {
            ofs << pos[i].x << " " << pos[i].y << " " << pos[i].z << " ";
            ofs << color[i].x << " " << color[i].y << " " << color[i].z << "\n";
        }

        ofs.close();
    }
}

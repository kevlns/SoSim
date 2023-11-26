//@author        : Long Shen
//@date          : 2023/11/21
//@description   :
//@version       : 1.0

#ifndef SOSIM_SHADER_HPP
#define SOSIM_SHADER_HPP

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

namespace SoSim::GUI {

    class Shader {
    public:
        unsigned int ID;

        // constructor generates the shader on the fly
        // ------------------------------------------------------------------------
        Shader(const char *vertexPath, const char *fragmentPath);

        // activate the shader
        // ------------------------------------------------------------------------
        inline void use();

        // utility uniform functions
        // ------------------------------------------------------------------------
        inline void setBool(const std::string &name, bool value) const;

        // ------------------------------------------------------------------------
        inline void setInt(const std::string &name, int value) const;

        // ------------------------------------------------------------------------
        void setFloat(const std::string &name, float value) const;

    private:
        // utility function for checking shader compilation/linking errors.
        // ------------------------------------------------------------------------
        void checkCompileErrors(unsigned int shader, std::string type);
    };
}

#endif // SOSIM_SHADER_HPP

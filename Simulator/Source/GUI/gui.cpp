//@author        : Long Shen
//@date          : 2023/11/21
//@description   :
//@version       : 1.0

#include "Public/GUI/gui.hpp"
#include "Private/GUI/widgets.hpp"
#include "Public/Shared/ModelUtils/model_builder.hpp"

namespace SoSim {

    Camera camera = Camera(glm::vec3(-5.f, 8.f, 35.0f), glm::vec3(0.0f, 35.0f, -8.f));
    const unsigned int SCR_WIDTH = 1024;
    const unsigned int SCR_HEIGHT = 720;

    static void glfw_error_callback(int error, const char *description) {
        fprintf(stderr, "GLFW Error %d: %s\n", error, description);
    }

    GUI::GUI() {
        if (!m_simulator)
            m_simulator = new Simulator;

//        if (!m_renderer)
//            m_renderer = new Renderer;
    }

    void GUI::initialize() {
        glfwSetErrorCallback(glfw_error_callback);
        if (!glfwInit())
            return;

        const char *glsl_version = "#version 130";
        glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
        glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);

        // Create window with graphics context
        m_main_window = glfwCreateWindow(1280, 720, "Dear ImGui GLFW+OpenGL3 example", nullptr, nullptr);
        if (m_main_window == nullptr)
            return;
        glfwMakeContextCurrent(m_main_window);
        glfwSwapInterval(1); // Enable vsync

        // Setup Dear ImGui context
        IMGUI_CHECKVERSION();
        ImGui::CreateContext();
        m_io = ImGui::GetIO();
        (void) m_io;
        m_io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard; // Enable Keyboard Controls
        m_io.ConfigFlags |= ImGuiConfigFlags_NavEnableGamepad; // Enable Gamepad Controls

        // Setup Dear ImGui style
        ImGui::StyleColorsDark();
        //ImGui::StyleColorsLight();

        // Setup Platform/Renderer backends
        ImGui_ImplGlfw_InitForOpenGL(m_main_window, true);
        ImGui_ImplOpenGL3_Init(glsl_version);

        m_io.Fonts->AddFontFromFileTTF("c:\\Windows\\Fonts\\STXINWEI.ttf", 15.0f);
    }

    void GUI::run(bool keepSilent) {
        m_keep_silent = keepSilent;
        if (keepSilent)
            runPure();
        else
            runGUI();
    }

    void GUI::terminate() {
        if (!m_keep_silent) {
            // Cleanup
            ImGui_ImplOpenGL3_Shutdown();
            ImGui_ImplGlfw_Shutdown();
            ImGui::DestroyContext();

            glfwDestroyWindow(m_main_window);
            glfwTerminate();
            m_main_window = nullptr;
        }

        m_simulator->terminate();
        delete m_simulator;
        m_simulator = nullptr;
    }

    void GUI::runPure() {
        // TODO
    }

    void GUI::runGUI() {
        initialize();

        if (!gladLoadGLLoader((GLADloadproc) glfwGetProcAddress)) {
            // 初始化失败的处理
            std::cerr << "Failed to initialize GLAD" << std::endl;
        }

        /**
         * test draw sphere
         */
        ObjectConfig objectConfig{};
        objectConfig.shape = "cube";
        objectConfig.lb = {-3, -3, -3};
        objectConfig.size = {6, 6, 6};
        objectConfig.color = {155, 155, 155};
        objectConfig.particle_radius = 0.05;

        auto pos = genParticleCube(&objectConfig);

        std::vector<float3> newPos;
        const double M_PI = 3.1415926;
        int latitudeBands = 5; // 纬度带数量
        int longitudeBands = 5; // 经度带数量
        float radius = 0.05; // 球体半径

        for (int latNumber = 0; latNumber <= latitudeBands; latNumber++) {
            float theta = latNumber * M_PI / latitudeBands;
            float sinTheta = sin(theta);
            float cosTheta = cos(theta);

            for (int longNumber = 0; longNumber <= longitudeBands; longNumber++) {
                float phi = longNumber * 2 * M_PI / longitudeBands;
                float sinPhi = sin(phi);
                float cosPhi = cos(phi);

                float x = cosPhi * sinTheta;
                float y = cosTheta;
                float z = sinPhi * sinTheta;

                // 将顶点坐标添加到顶点数组
                newPos.push_back({x, y, z});
            }
        }

        std::vector<unsigned int> indices;
        for (int latNumber = 0; latNumber < latitudeBands; latNumber++) {
            for (int longNumber = 0; longNumber < longitudeBands; longNumber++) {
                unsigned int first = (latNumber * (longitudeBands + 1)) + longNumber;
                unsigned int second = first + longitudeBands + 1;

                indices.push_back(first);
                indices.push_back(second);
                indices.push_back(first + 1);

                indices.push_back(second);
                indices.push_back(second + 1);
                indices.push_back(first + 1);
            }
        }

        unsigned int VBO, VAO, EBO;
        glGenVertexArrays(1, &VAO);
        glGenBuffers(1, &VBO);
        glGenBuffers(1, &EBO);

        glBindVertexArray(VAO);

        glBindBuffer(GL_ARRAY_BUFFER, VBO);
        glBufferData(GL_ARRAY_BUFFER, newPos.size() * sizeof(float), newPos.data(), GL_STATIC_DRAW);

        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(unsigned int), indices.data(), GL_STATIC_DRAW);

// 设置顶点属性指针
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void *) 0);
        glEnableVertexAttribArray(0);

        glBindBuffer(GL_ARRAY_BUFFER, 0);
        glBindVertexArray(0);

//        unsigned int VAO1, VBO1;
//        glGenVertexArrays(1, &VAO1);
//        glGenBuffers(1, &VBO1);
//
//        glBindVertexArray(VAO1);
//
//        glBindBuffer(GL_ARRAY_BUFFER, VBO1);
//        auto pointNum = pos.size();
//        glBufferData(GL_ARRAY_BUFFER, pointNum * sizeof(float3), pos.data(), GL_STATIC_DRAW);
//
//        // position attribute
//        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void *) 0);
//        glEnableVertexAttribArray(0);
//
//        glBindBuffer(GL_ARRAY_BUFFER, 0);

        glm::mat4 model = glm::mat4(1.0f);
        glm::mat4 projection = glm::perspective(glm::radians(camera.Zoom), (float) SCR_WIDTH / (float) SCR_HEIGHT,
                                                0.1f,
                                                100.0f);

        Shader shader("test_1.vert", "test_1.frag");
        shader.use();
        shader.setMat4("model", model);
        shader.setMat4("projection", projection);
        shader.setVec3("objectColor", {0.2, 0.3, 0.4});
        shader.setVec3("lightColor", {1, 1, 1});
        shader.setVec3("lightPos", {5, 5, 5});
        shader.setVec3("viewPos", camera.Position);


        // Our state
        bool show_demo_window = true;
        bool show_another_window = false;
        ImVec4 clear_color = ImVec4(0.2f, 0.2f, 0.2f, 1.00f);

        glEnable(GL_POINT_SPRITE);
        // Main loop
        while (!glfwWindowShouldClose(m_main_window)) {
            // Poll and handle events (inputs, window resize, etc.)
            // You can read the io.WantCaptureMouse, io.WantCaptureKeyboard flags to tell if dear imgui wants to use your inputs.
            // - When io.WantCaptureMouse is true, do not dispatch mouse input data to your main application, or clear/overwrite your copy of the mouse data.
            // - When io.WantCaptureKeyboard is true, do not dispatch keyboard input data to your main application, or clear/overwrite your copy of the keyboard data.
            // Generally you may always pass all inputs to dear imgui, and hide them from your application based on those two flags.
            glfwPollEvents();

            // Start the Dear ImGui frame
            ImGui_ImplOpenGL3_NewFrame();
            ImGui_ImplGlfw_NewFrame();
            ImGui::NewFrame();


            ShowMainWindow(m_simulator);

//            if (show_demo_window)
//                ImGui::ShowDemoWindow(&show_demo_window);

            // Rendering
            ImGui::Render();
            int display_w, display_h;
            glfwGetFramebufferSize(m_main_window, &display_w, &display_h);
            glViewport(0, 0, display_w, display_h);
            glClearColor(clear_color.x * clear_color.w, clear_color.y * clear_color.w, clear_color.z * clear_color.w,
                         clear_color.w);
            glClear(GL_COLOR_BUFFER_BIT);
            ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

            glm::mat4 view = camera.GetViewMatrix();
            shader.use();
            shader.setMat4("view", view);
            shader.setVec4("pureColor", {0.8, 0.3, 0.1, 1.0});

            glBindVertexArray(VAO);
            glDrawElements(GL_TRIANGLES, indices.size(), GL_UNSIGNED_INT, 0);
            glBindVertexArray(0);
//            glBindBuffer(GL_ARRAY_BUFFER, VBO1);
//            glDrawArrays(GL_POINTS, 0, pointNum);
//            glBindBuffer(GL_ARRAY_BUFFER, 0);

            glfwSwapBuffers(m_main_window);
        }

        glDeleteVertexArrays(1, &VAO);
        glDeleteBuffers(1, &VBO);
        glDeleteBuffers(1, &EBO);
    }
}

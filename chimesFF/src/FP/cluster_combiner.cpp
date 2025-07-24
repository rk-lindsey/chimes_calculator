#include <iostream>
#include <filesystem>
#include <regex>
#include <map>
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>

namespace fs = std::filesystem;

int main() {
    const std::vector<int> target_bodies = {2, 3, 4};

    for (int body : target_bodies) {
        std::regex filename_pattern(R"(^(\d+)\.(\d+)\.)" + std::to_string(body) + R"(b_clusters\.txt$)");
        std::map<int, std::vector<std::pair<int, fs::path>>> frame_files;

        // Collect files
        for (const auto& entry : fs::directory_iterator(".")) {
            std::string filename = entry.path().filename().string();
            std::smatch match;
            
            if (std::regex_match(filename, match, filename_pattern)) {
                try {
                    int frame = std::stoi(match[1]);
                    int rank = std::stoi(match[2]);
                    frame_files[frame].emplace_back(rank, entry.path());
                } catch (...) {}
            }
        }

        // Process frames
        for (const auto& [frame, files] : frame_files) {
            auto sorted_files = files;
            std::sort(sorted_files.begin(), sorted_files.end(),
                [](const auto& a, const auto& b) { return a.first < b.first; });

            std::ostringstream oss;
            oss << frame << "." << body << "b_combined.txt";
            fs::path output_path = oss.str();

            // 1. Clear any existing output stream errors
            std::ofstream output;
            output.clear();
            
            // 2. Open in append mode to preserve existing content
            output.open(output_path, std::ios::binary | std::ios::app);
            if (!output.is_open()) {
                std::cerr << "CRITICAL ERROR: Failed to create " << output_path << "\n";
                continue;
            }

            size_t initial_position = output.tellp();
            size_t total_bytes = 0;

            for (const auto& [rank, path] : sorted_files) {
                std::error_code ec;
                size_t file_size = fs::file_size(path, ec);
                
                std::cout << "Processing: " << path.filename() 
                          << " (" << file_size << " bytes)\n";

                std::ifstream input(path, std::ios::binary);
                if (!input) {
                    std::cerr << "  FILE ERROR: Cannot open input file\n";
                    continue;
                }

                // 3. Read file contents into memory first
                std::string content((std::istreambuf_iterator<char>(input)), 
                                     std::istreambuf_iterator<char>());
                
                // 4. Verify we actually got the data
                if (content.size() != file_size) {
                    std::cerr << "  READ ERROR: Only read " << content.size() 
                              << "/" << file_size << " bytes\n";
                    continue;
                }

                // 5. Write from memory buffer
                output.write(content.data(), content.size());
                
                // 6. Immediate flush and check
                output.flush();
                if (!output) {
                    std::cerr << "  WRITE ERROR: Failed to write contents\n";
                    output.clear();
                    continue;
                }

                total_bytes += content.size();
            }

            // 7. Final verification
            size_t final_size = static_cast<size_t>(output.tellp()) - initial_position;
            std::cout << "VERIFICATION: " << output_path 
                      << " should be " << total_bytes << " bytes, "
                      << "actually " << final_size << " bytes\n\n";

            output.close();
        }
    }

    return 0;
}
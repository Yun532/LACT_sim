#include "core/Sha256.hpp"

#include <filesystem>
#include <fstream>
#include <iostream>

int main()
{
    const std::string expected_empty =
        "e3b0c44298fc1c149afbf4c8996fb924"
        "27ae41e4649b934ca495991b7852b855";
    const std::string expected_abc =
        "ba7816bf8f01cfea414140de5dae2223"
        "b00361a396177a9cb410ff61f20015ad";
    if (lact::sha256String("") != expected_empty ||
        lact::sha256String("abc") != expected_abc) {
        std::cerr << "SHA-256 known-answer test failed\n";
        return 1;
    }

    const auto path =
        std::filesystem::temp_directory_path() / "lact_sha256_abc.txt";
    {
        std::ofstream output(path, std::ios::binary);
        output << "abc";
    }
    const auto file_hash = lact::sha256File(path.string());
    std::filesystem::remove(path);
    if (file_hash != expected_abc) {
        std::cerr << "SHA-256 file hashing test failed\n";
        return 1;
    }
    return 0;
}

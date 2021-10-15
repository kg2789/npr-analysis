#include <algorithm>
#include <string>

#include <tclap/CmdLine.h>

#include "format.h"
#include "io.h"
#include "../data.h"

int main(int argc, char **argv) {
    Format input_format, output_format;
    std::string input_dir, output_dir;

    // Using tclap to parse command-line arguments; see here for more detail:
    // https://stackoverflow.com/questions/865668/how-to-parse-command-line-arguments-in-c
    try {
        TCLAP::CmdLine cmd("Data conversion program for NPR data");

        TCLAP::ValueArg<std::string>
            inputArg("i", "input", "Input directory",
                    true, "", "filename", cmd);
        TCLAP::ValueArg<std::string>
            outputArg("o", "output", "Output directory",
                    true, "", "filename", cmd);
        TCLAP::ValueArg<std::string>
            inputFormatArg("f", "input-format", "Input format",
                    false, "raw_v2", "format", cmd);
        TCLAP::ValueArg<std::string>
            outputFormatArg("r", "output-format", "Output format",
                    false, "rabbott", "format", cmd);

        cmd.parse(argc, argv);

        input_dir = inputArg.getValue();
        output_dir = outputArg.getValue();
        input_format  = parse_format(inputFormatArg.getValue());
        output_format = parse_format(outputFormatArg.getValue());
    }
    catch (TCLAP::ArgException &e) { // catch any exceptions
        std::cerr << "error: " << e.error()
            << " for arg " << e.argId() << std::endl;
    }

    std::cout << "Reading data from '" << input_dir << "'"
        << "in format " << input_format << std::endl;
    RunData data = read_data(input_dir, input_format, false);
    std::cout << "Writing data to '" << output_dir << "'"
        << "in format " << output_format << std::endl;
    write_data(output_dir, data, output_format);
    /*read_trajectory(input_dir, data, 1500, input_format);
    std::cout << "Writing tranformed data" << std::endl;
    write_trajectory(output_dir, data, 1500, output_format);

    std::cout << "Data norm: " << l2_norm(*data) << std::endl;

    std::cout << "Re-reading written data to test fidelity" << std::endl;
    TrajectoryData* read_data = new TrajectoryData();
    read_trajectory(output_dir, read_data, 1500, output_format);

    FourQOpFourQExtMeasurements& raw_measurments = data->fourq_op.fourq_ext;
    FourQOpFourQExtMeasurements& written_measurments = read_data->fourq_op.fourq_ext;

    for (int a = 0; a < NUM_LOOP_STRUCTURES; a++) {
        for (int b = 0; b < NUM_COLOR_STRUCTURES; b++) {
            for (int c = 0; c < NUM_SPIN_STRUCTURES; c++) {
                for (int d = 0; d < NUM_SPIN_STRUCTURES; d++) {
                    assert(raw_measurments.results[a][b][c][d]
                            == written_measurments.results[a][b][c][d]);
                }
            }
        }
    }
    delete data;
    delete read_data;*/
}

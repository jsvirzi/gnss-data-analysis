#include<cstdio>
#include<cstdlib>
#include <iostream>
#include <string>

#include <iostream>
#include <vector>
#include <sys/stat.h>
#include <fcntl.h>

#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <sstream>
#include <fstream>
#include <inttypes.h>
#include <stdint.h>
#include <cinttypes>
#include <iomanip>

#include <TApplication.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TMarker.h>
#include <TMath.h>
#include <TFile.h>
#include <TH1D.h>

std::vector<std::string> splitFields(const std::string &input_string, const std::string &delimiters, int option = 0) {
    std::vector<std::string> fields;
    std::string empty_string;
    std::string str = input_string;
    size_t pos = (option == 0) ? str.find_first_of(delimiters) : str.find(delimiters);
    size_t add_length = (option == 0) ? 1 : delimiters.size();
    while (pos != std::string::npos) {
        if (pos > 0) {
            fields.push_back(str.substr(0, pos));
        } else if (pos == 0) {
            fields.push_back(empty_string);
        }
        str = str.substr(pos + add_length);
        pos = (option == 0) ? str.find_first_of(delimiters) : str.find(delimiters);
    }
    if (str.length()) {
        fields.push_back(str);
    } else {
        fields.push_back(empty_string);
    }
    return fields;
}

enum {
    EgoTimestampCsvIndex = 0,
    EgoTimestampLatIndex = 1,
    EgoTimestampLonIndex = 2,
    EgoTimestampAltIndex = 3,
    EgoTimestampLatSigmaIndex = 4,
    EgoTimestampLonSigmaIndex = 5,
    EgoTimestampAltSigmaIndex = 6,
};

enum {
    IpaTimestampCsvIndex = 0,
    IpaTimestampLatIndex = 2,
    IpaTimestampLonIndex = 3,
    IpaTimestampAltIndex = 4,
    IpaTimestampLatSigmaIndex = 5,
    IpaTimestampLonSigmaIndex = 6,
    IpaTimestampAltSigmaIndex = 7,
};

#if 0
(0.0001,0.0001) = (8.824m, 11.112m)
#endif
constexpr double kLatCalibration = 88240.0;
constexpr double kLonCalibration = 111120.0;

int main(int argc, char **argv) {
    uint64_t timestamp;
    double lat, lat_sigma, lon, lon_sigma, alt, alt_sigma;
    std::ifstream *ifs;
    std::vector<double> *v_lat, *v_lon, *v_alt, *v_lat_sigma, *v_lon_sigma, *v_alt_sigma;

    TApplication the_app("analysis", 0, 0);

    /*** read in data from out device ***/
    std::string iGpsFileIpl = "/home/ubuntu/data/analysis/gps_parsed.csv";
    std::ifstream ifstream_ego;
    ifstream_ego.open(iGpsFileIpl);
    ifs = &ifstream_ego;
    if (ifs->is_open() == false) {
        printf("boo hoo\n");
        return 1;
    }

    std::vector<double> v_ego_lat, v_ego_lon, v_ego_alt, v_ego_lat_sigma, v_ego_lon_sigma, v_ego_alt_sigma;
    v_lat = &v_ego_lat; v_lat_sigma = &v_ego_lat_sigma;
    v_lon = &v_ego_lon; v_lon_sigma = &v_ego_lon_sigma;
    v_alt = &v_ego_alt; v_alt_sigma = &v_ego_alt_sigma;
    int n_gps_ego = 0;
    for (std::string line; std::getline(*ifs, line);) {
        std::vector<std::string> fields = splitFields(std::string(line), ",");
        sscanf(fields[EgoTimestampCsvIndex].c_str(), "%lu", &timestamp);
        sscanf(fields[EgoTimestampLatIndex].c_str(), "%lf", &lat);
        sscanf(fields[EgoTimestampLonIndex].c_str(), "%lf", &lon);
        sscanf(fields[EgoTimestampAltIndex].c_str(), "%lf", &alt);
        sscanf(fields[EgoTimestampLatSigmaIndex].c_str(), "%lf", &lat_sigma);
        sscanf(fields[EgoTimestampLonSigmaIndex].c_str(), "%lf", &lon_sigma);
        sscanf(fields[EgoTimestampAltSigmaIndex].c_str(), "%lf", &alt_sigma);
        v_lat->push_back(lat);
        v_lat_sigma->push_back(lat_sigma);
//        printf("line = %lu %lf +/- %lf %lf +/- %lf\n", timestamp, lat, lat_sigma, lon, lon_sigma);
        ++n_gps_ego;
    }

    printf("n = %d\n", n_gps_ego);

    uint64_t *ego_timestamp = new uint64_t [ n_gps_ego ];
    double *ego_lat = new double [ n_gps_ego ];
    double *ego_lon = new double [ n_gps_ego ];
    double *ego_alt = new double [ n_gps_ego ];
    double *ego_lat_sigma = new double [ n_gps_ego ];
    double *ego_lon_sigma = new double [ n_gps_ego ];
    double *ego_alt_sigma = new double [ n_gps_ego ];

    ifstream_ego.clear();
    ifstream_ego.seekg(0, std::ios::beg); /* rewind */

    int i_gps = 0;
    for (std::string line; std::getline(ifstream_ego, line);++i_gps) {
        std::vector<std::string> fields = splitFields(std::string(line), ",");
        sscanf(fields[EgoTimestampCsvIndex].c_str(), "%lu", &timestamp);
        sscanf(fields[EgoTimestampLatIndex].c_str(), "%lf", &lat);
        sscanf(fields[EgoTimestampLonIndex].c_str(), "%lf", &lon);
        sscanf(fields[EgoTimestampAltIndex].c_str(), "%lf", &alt);
        sscanf(fields[EgoTimestampLatSigmaIndex].c_str(), "%lf", &lat_sigma);
        sscanf(fields[EgoTimestampLonSigmaIndex].c_str(), "%lf", &lon_sigma);
        sscanf(fields[EgoTimestampAltSigmaIndex].c_str(), "%lf", &alt_sigma);
        ego_timestamp[i_gps] = timestamp;
        ego_lat[i_gps] = lat * kLatCalibration;
        ego_lon[i_gps] = lon * kLonCalibration;
        ego_alt[i_gps] = alt;
        ego_lat_sigma[i_gps] = lat_sigma;
        ego_lon_sigma[i_gps] = lon_sigma;
        ego_alt_sigma[i_gps] = alt_sigma;
        printf("line = %lu %lf +/- %lf %lf +/- %lf\n", timestamp, lat, lat_sigma, lon, lon_sigma);
    }

    ifstream_ego.close();

    /*** read in IPA data ***/
    std::string iGpsFileIpa = "/home/ubuntu/data/analysis/ipa_ipl_parsed.csv";
    std::ifstream ifstream_ipa;
    ifstream_ipa.open(iGpsFileIpa);
    if (ifstream_ipa.is_open() == false) {
        printf("boo hoo\n");
    }

    int n_gps_ipa = 0;

    for (std::string line; std::getline(ifstream_ipa, line);) {
        std::vector<std::string> fields = splitFields(std::string(line), ",");
        sscanf(fields[IpaTimestampCsvIndex].c_str(), "%lu", &timestamp);
        sscanf(fields[IpaTimestampLatIndex].c_str(), "%lf", &lat);
        sscanf(fields[IpaTimestampLonIndex].c_str(), "%lf", &lon);
        sscanf(fields[IpaTimestampAltIndex].c_str(), "%lf", &alt);
        sscanf(fields[IpaTimestampLatSigmaIndex].c_str(), "%lf", &lat_sigma);
        sscanf(fields[IpaTimestampLonSigmaIndex].c_str(), "%lf", &lon_sigma);
        sscanf(fields[IpaTimestampAltSigmaIndex].c_str(), "%lf", &alt_sigma);
//        printf("line = %lu %lf +/- %lf %lf +/- %lf\n", timestamp, lat, lat_sigma, lon, lon_sigma);
        ++n_gps_ipa;
    }

    printf("n = %d\n", n_gps_ipa);

    uint64_t *ipa_timestamp = new uint64_t [ n_gps_ipa ];
    double *ipa_lat = new double [ n_gps_ipa ];
    double *ipa_lon = new double [ n_gps_ipa ];
    double *ipa_alt = new double [ n_gps_ipa ];
    double *ipa_lat_sigma = new double [ n_gps_ipa ];
    double *ipa_lon_sigma = new double [ n_gps_ipa ];
    double *ipa_alt_sigma = new double [ n_gps_ipa ];

    ifstream_ipa.clear();
    ifstream_ipa.seekg(0, std::ios::beg); /* rewind */

    i_gps = 0;
    for (std::string line; std::getline(ifstream_ipa, line);) {
        std::vector<std::string> fields = splitFields(std::string(line), ",");
        sscanf(fields[IpaTimestampCsvIndex].c_str(), "%lu", &timestamp);
        sscanf(fields[IpaTimestampLatIndex].c_str(), "%lf", &lat);
        sscanf(fields[IpaTimestampLonIndex].c_str(), "%lf", &lon);
        sscanf(fields[IpaTimestampAltIndex].c_str(), "%lf", &alt);
        sscanf(fields[IpaTimestampLatSigmaIndex].c_str(), "%lf", &lat_sigma);
        sscanf(fields[IpaTimestampLonSigmaIndex].c_str(), "%lf", &lon_sigma);
        sscanf(fields[IpaTimestampAltSigmaIndex].c_str(), "%lf", &alt_sigma);
        ipa_timestamp[i_gps] = timestamp;
        ipa_lat[i_gps] = lat * kLatCalibration;
        ipa_lon[i_gps] = lon * kLonCalibration;
        ipa_alt[i_gps] = alt;
        ipa_lat_sigma[i_gps] = lat_sigma;
        ipa_lon_sigma[i_gps] = lon_sigma;
        ipa_alt_sigma[i_gps] = alt_sigma;
        printf("line = %lu %lf +/- %lf %lf +/- %lf\n", timestamp, lat, lat_sigma, lon, lon_sigma);
        ++i_gps;
    }

    ifstream_ipa.close();

    /* remove offsets */
    double ego_lat_mean = 0.0, ego_lon_mean = 0.0, ego_alt_mean = 0.0;
    for (i_gps = 0; i_gps < n_gps_ego; ++i_gps) {
        ego_lat_mean += ego_lat[i_gps];
        ego_lon_mean += ego_lon[i_gps];
        ego_alt_mean += ego_alt[i_gps];
    }
    ego_lat_mean /= n_gps_ego;
    ego_lon_mean /= n_gps_ego;
    ego_alt_mean /= n_gps_ego;
    for (i_gps = 0; i_gps < n_gps_ego; ++i_gps) {
        ego_lat[i_gps] -= ego_lat_mean;
        ego_lon[i_gps] -= ego_lon_mean;
        ego_alt[i_gps] -= ego_alt_mean;
    }

    double *time_ego = new double [ n_gps_ego ];
    for (i_gps = 0; i_gps < n_gps_ego; ++i_gps) {
        uint64_t dtime = ego_timestamp[i_gps] - ego_timestamp[0];
        time_ego[i_gps] = (double) dtime / 1000000.f;
    }

    /* remove offsets */
    double ipa_lat_mean = 0.0, ipa_lon_mean = 0.0, ipa_alt_mean = 0.0;
    for (i_gps = 0; i_gps < n_gps_ipa; ++i_gps) {
        ipa_lat_mean += ipa_lat[i_gps];
        ipa_lon_mean += ipa_lon[i_gps];
        ipa_alt_mean += ipa_alt[i_gps];
    }
    ipa_lat_mean /= n_gps_ipa;
    ipa_lon_mean /= n_gps_ipa;
    ipa_alt_mean /= n_gps_ipa;
    for (i_gps = 0; i_gps < n_gps_ipa; ++i_gps) {
        ipa_lat[i_gps] -= ipa_lat_mean;
        ipa_lon[i_gps] -= ipa_lon_mean;
        ipa_alt[i_gps] -= ipa_alt_mean;
    }

    double *time_ipa = new double [ n_gps_ipa ];
    for (i_gps = 0; i_gps < n_gps_ipa; ++i_gps) {
        uint64_t dtime = ipa_timestamp[i_gps] - ipa_timestamp[0];
        time_ipa[i_gps] = (double) dtime / 1000.f;
    }

    double *dlat = new double [ n_gps_ego ];
    double *ddlat = new double [ n_gps_ego ];
    double *ddt = new double [ n_gps_ego ];
    for (i_gps = 0; i_gps < n_gps_ego; ++i_gps) {
        dlat[i_gps] = ego_lat[i_gps] - ipa_lat[i_gps];
        double elat_ego = ego_lat_sigma[i_gps];
        double elat_ipa = ipa_lat_sigma[i_gps];
        ddlat[i_gps] = TMath::Sqrt(elat_ego * elat_ego + elat_ipa * elat_ipa);
        ddt[i_gps] = 0.0;
    }

    /*** get min/max ***/
    double min_lat = 99999.0, max_lat = -999999.0;
    for (i_gps = 0; i_gps < n_gps_ego; ++i_gps) {
        double lat = ego_lat[i_gps];
        printf("ego: %d time=%lf lat=%lf\n", i_gps, time_ego[i_gps], lat);
        if (lat > max_lat) { max_lat = lat; }
        if (lat < min_lat) { min_lat = lat; }
    }

    for (i_gps = 0; i_gps < n_gps_ipa; ++i_gps) {
        double lat = ipa_lat[i_gps];
        printf("ipa: %d time=%lf lat=%lf\n", i_gps, time_ipa[i_gps], lat);
        if (lat > max_lat) { max_lat = lat; }
        if (lat < min_lat) { min_lat = lat; }
    }

    double min_dlat = 99999.0, max_dlat = -999999.0;
    for (i_gps = 0; i_gps < n_gps_ego; ++i_gps) {
        double lat = dlat[i_gps];
        if (lat > max_dlat) { max_dlat = lat; }
        if (lat < min_dlat) { min_dlat = lat; }
    }

    printf("lat min/max = %lf/%lf\n", min_lat, max_lat);

    int n_graph = n_gps_ego;
    n_graph = 200;
    TGraphErrors *graph_dlat = new TGraphErrors(n_graph, time_ego, dlat, ddt, ddlat);
    TGraphErrors *graph_lat_ego = new TGraphErrors(n_gps_ego, time_ego, ego_lat, ddt, ego_lat_sigma);
    TGraphErrors *graph_lat_ipa = new TGraphErrors(n_gps_ipa, time_ipa, ipa_lat, ddt, ipa_lat_sigma);
    TH1D *h_dlat = new TH1D("dlat", "dlat", n_graph, time_ego[0], time_ego[n_graph]);

    TFile fp("/home/ubuntu/data/analysis/analysis.root", "recreate");
    graph_dlat->Write("dlat");
    graph_lat_ego->Write("lat_ego");
    graph_lat_ipa->Write("lat_ipa");
    fp.Write();
    fp.Close();

    TCanvas canvas("main", "main", 1200, 800);
    TGraphErrors *graph;

    graph = graph_dlat;
    graph->SetMinimum(min_dlat);
    graph->SetMaximum(max_dlat);
    graph->SetMarkerStyle(kFullCircle);
    graph->SetFillColor(kYellow);
    graph->SetFillStyle(kSolid);
    graph->Draw("AP E0");
//    graph = graph_lat_ego;
//    graph->SetMarkerStyle(kFullTriangleUp);
//    graph->Draw("P");
    canvas.Draw();
    canvas.Update();
    canvas.WaitPrimitive();

    getchar();

    return 0;
}

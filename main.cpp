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
#include <TLegend.h>

void runExperiment();

#if 0
double rms(double *x, int n) {
    double acc0 = 0.0, acc1 = 0.0, acc2 = 0.0;
    for (int i=0;i<n;++i) {
        double a = *x++;
        acc0 += 1.0;
        acc1 += a;
        acc2 += a * a;
    }
    double mean = acc1 / acc0;
    double std_dev = sqrt(fabs(acc2 / acc0 - mean * mean));
    return std_dev;
}
#else
double rms(double *x, int n) {
    double acc0 = 0.0, acc1 = 0.0, acc2 = 0.0;
    for (int i=1;i<n;++i) {
        double a = x[i] - x[i-1];
        acc0 += 1.0;
        acc1 += a;
        acc2 += a * a;
    }
    return acc2 / acc0;
}
#endif

void purdifyPlots(TH1D *h1, TH1D *h2, const char *title, const char *x_title) {
    TH1D *h;
    h = h1;
    h->SetFillColor(kYellow);
    h->SetLineColor(kBlack);
    h->SetLineWidth(3.0);
    h->SetStats(kFALSE);
    h->GetXaxis()->SetTitle(x_title);
    h->GetYaxis()->SetTitle("arbitrary units");
    h->GetXaxis()->CenterTitle(kTRUE);
    h->GetYaxis()->CenterTitle(kTRUE);
    h->SetTitle(title);

    h = h2;
    h->SetLineColor(kRed);
    h->SetLineWidth(3.0);
    h->SetMarkerStyle(kFullCircle);
    h->SetMarkerColor(kRed);
    h->SetMarkerSize(1.0);
    h->SetStats(kFALSE);
    h->GetXaxis()->SetTitle(x_title);
    h->GetYaxis()->SetTitle("arbitrary units");
    h->GetXaxis()->CenterTitle(kTRUE);
    h->GetYaxis()->CenterTitle(kTRUE);
    h->SetTitle(title);
}

std::string data_directory = "/home/jsvirzi/projects/gnss-data-analysis/data";
// std::string data_directory = "/home/ubuntu/projects/gnss-data-analysis/data";

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
    std::vector<uint64_t> *v_timestamp;
    int i_gps;
    TLegend *legend;

    TApplication the_app("analysis", 0, 0);

//    runExperiment();

    /*** read in data from out device ***/
    std::string iGpsFileIpl = data_directory + "/gps_parsed.csv";
    std::ifstream ifstream_ego;
    ifstream_ego.open(iGpsFileIpl);
    ifs = &ifstream_ego;
    if (ifs->is_open() == false) {
        printf("boo hoo\n");
        return 1;
    }

    std::vector<double> v_ego_lat, v_ego_lon, v_ego_alt, v_ego_lat_sigma, v_ego_lon_sigma, v_ego_alt_sigma;
    std::vector<uint64_t> v_ego_timestamp;
    v_lat = &v_ego_lat; v_lat_sigma = &v_ego_lat_sigma;
    v_lon = &v_ego_lon; v_lon_sigma = &v_ego_lon_sigma;
    v_alt = &v_ego_alt; v_alt_sigma = &v_ego_alt_sigma;
    v_timestamp = &v_ego_timestamp;
    for (std::string line; std::getline(*ifs, line);) {
        std::vector<std::string> fields = splitFields(std::string(line), ",");
        sscanf(fields[EgoTimestampCsvIndex].c_str(), "%lu", &timestamp);
        sscanf(fields[EgoTimestampLatIndex].c_str(), "%lf", &lat);
        sscanf(fields[EgoTimestampLonIndex].c_str(), "%lf", &lon);
        sscanf(fields[EgoTimestampAltIndex].c_str(), "%lf", &alt);
        sscanf(fields[EgoTimestampLatSigmaIndex].c_str(), "%lf", &lat_sigma);
        sscanf(fields[EgoTimestampLonSigmaIndex].c_str(), "%lf", &lon_sigma);
        sscanf(fields[EgoTimestampAltSigmaIndex].c_str(), "%lf", &alt_sigma);
        lat *= kLatCalibration; lon *= kLonCalibration;
        v_timestamp->push_back(timestamp);
        v_lat->push_back(lat); v_lat_sigma->push_back(lat_sigma);
        v_lon->push_back(lon); v_lon_sigma->push_back(lon_sigma);
        v_alt->push_back(alt); v_alt_sigma->push_back(alt_sigma);
//        printf("line = %lu %lf +/- %lf %lf +/- %lf\n", timestamp, lat, lat_sigma, lon, lon_sigma);
    }

    int n_gps_ego = v_timestamp->size();
    printf("ego: n = %d\n", n_gps_ego);

    uint64_t *ego_timestamp = v_ego_timestamp.data();
    double *ego_lat = v_ego_lat.data();
    double *ego_lon = v_ego_lon.data();
    double *ego_alt = v_ego_alt.data();
    double *ego_lat_sigma = v_ego_lat_sigma.data();
    double *ego_lon_sigma = v_ego_lon_sigma.data();
    double *ego_alt_sigma = v_ego_alt_sigma.data();

    ifs->clear(); /* we read to the end, so need to clear EOF bit and then seekg() should work */
    ifs->seekg(0, std::ios::beg); /* rewind */
    ifs->close();

    /*** read in IPA data ***/
    std::string iGpsFileIpa = data_directory + "/ipa_ipl_parsed.csv";
    std::ifstream ifstream_ipa;
    ifstream_ipa.open(iGpsFileIpa);
    ifs = &ifstream_ipa;
    if (ifs->is_open() == false) {
        printf("boo hoo\n");
        return 1;
    }

    std::vector<double> v_ipa_lat, v_ipa_lon, v_ipa_alt, v_ipa_lat_sigma, v_ipa_lon_sigma, v_ipa_alt_sigma;
    std::vector<uint64_t> v_ipa_timestamp;
    v_lat = &v_ipa_lat; v_lat_sigma = &v_ipa_lat_sigma;
    v_lon = &v_ipa_lon; v_lon_sigma = &v_ipa_lon_sigma;
    v_alt = &v_ipa_alt; v_alt_sigma = &v_ipa_alt_sigma;
    v_timestamp = &v_ipa_timestamp;
    for (std::string line; std::getline(*ifs, line);) {
        std::vector<std::string> fields = splitFields(std::string(line), ",");
        sscanf(fields[IpaTimestampCsvIndex].c_str(), "%lu", &timestamp);
        sscanf(fields[IpaTimestampLatIndex].c_str(), "%lf", &lat);
        sscanf(fields[IpaTimestampLonIndex].c_str(), "%lf", &lon);
        sscanf(fields[IpaTimestampAltIndex].c_str(), "%lf", &alt);
        sscanf(fields[IpaTimestampLatSigmaIndex].c_str(), "%lf", &lat_sigma);
        sscanf(fields[IpaTimestampLonSigmaIndex].c_str(), "%lf", &lon_sigma);
        sscanf(fields[IpaTimestampAltSigmaIndex].c_str(), "%lf", &alt_sigma);
        lat *= kLatCalibration; lon *= kLonCalibration;
        v_timestamp->push_back(timestamp);
        v_lat->push_back(lat); v_lat_sigma->push_back(lat_sigma);
        v_lon->push_back(lon); v_lon_sigma->push_back(lon_sigma);
        v_alt->push_back(alt); v_alt_sigma->push_back(alt_sigma);
//        printf("line = %lu %lf +/- %lf %lf +/- %lf\n", timestamp, lat, lat_sigma, lon, lon_sigma);
    }

    int n_gps_ipa = v_timestamp->size();

    printf("ipa: n = %d\n", n_gps_ipa);

    uint64_t *ipa_timestamp = v_ipa_timestamp.data();
    double *ipa_lat = v_ipa_lat.data();
    double *ipa_lon = v_ipa_lon.data();
    double *ipa_alt = v_ipa_alt.data();
    double *ipa_lat_sigma = v_ipa_lat_sigma.data();
    double *ipa_lon_sigma = v_ipa_lon_sigma.data();
    double *ipa_alt_sigma = v_ipa_alt_sigma.data();

    ifs->clear(); /* we read to the end, so need to clear EOF bit and then seekg() should work */
    ifs->seekg(0, std::ios::beg); /* rewind */
    ifs->close();

    /* remove offsets from ego sample */
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

    /* remove offsets from ipa sample (not identical format to ego sample) */
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

    int n_gps = (n_gps_ego < n_gps_ipa) ? n_gps_ego : n_gps_ipa;
    double *dlat = new double [ n_gps ]; /* diff between ego and ipa */
    double *dlat_rms = new double [ n_gps ]; /* rms diff between ego and ipa */
    double *ddt = new double [ n_gps ]; /* time errors */
    for (i_gps = 0; i_gps < n_gps; ++i_gps) {
        dlat[i_gps] = ego_lat[i_gps] - ipa_lat[i_gps];
        double elat_ego = ego_lat_sigma[i_gps];
        double elat_ipa = ipa_lat_sigma[i_gps];
        dlat_rms[i_gps] = TMath::Sqrt(elat_ego * elat_ego + elat_ipa * elat_ipa);
        ddt[i_gps] = 0.0; /* zero for now */
    }

    /*** get min/max across both samples ***/
    double min_lat = 99999.0, max_lat = -999999.0;
    double max_ego_lat_sigma = -999999.0, max_ipa_lat_sigma = -999999.0;
    for (i_gps = 0; i_gps < n_gps_ego; ++i_gps) {
        double lat = ego_lat[i_gps];
//        printf("ego: %d time=%lf lat=%lf\n", i_gps, time_ego[i_gps], lat);
        if (lat > max_lat) { max_lat = lat; }
        if (lat < min_lat) { min_lat = lat; }
        double sigma = ego_lat_sigma[i_gps];
        if (sigma > max_ego_lat_sigma) { max_ego_lat_sigma = sigma; }
    }

    for (i_gps = 0; i_gps < n_gps_ipa; ++i_gps) {
        double lat = ipa_lat[i_gps];
//        printf("ipa: %d time=%lf lat=%lf\n", i_gps, time_ipa[i_gps], lat);
        if (lat > max_lat) { max_lat = lat; }
        if (lat < min_lat) { min_lat = lat; }
        double sigma = ipa_lat_sigma[i_gps];
        if (sigma > max_ipa_lat_sigma) { max_ipa_lat_sigma = sigma; }
    }

    double rms_ipa = rms(ipa_lat, n_gps_ipa);

    double min_dlat = 99999.0, max_dlat = -999999.0;
    for (i_gps = 0; i_gps < n_gps_ego; ++i_gps) {
        double lat = dlat[i_gps];
        if (lat > max_dlat) { max_dlat = lat; }
        if (lat < min_dlat) { min_dlat = lat; }
    }

    double rms_ego = rms(ego_lat, n_gps_ego);

    printf("lat min/max = %lf/%lf\n", min_lat, max_lat);
    printf("rms %lf(ego) %lf(ipa)\n", rms_ego, rms_ipa);

    /* graphs */
    int n_graph = n_gps_ego;
    n_graph = 200;
    TGraphErrors *graph_dlat = new TGraphErrors(n_graph, time_ego, dlat, ddt, dlat_rms);
    n_graph = n_gps_ego;
//    n_graph = 200;
    TGraphErrors *graph_lat_ego = new TGraphErrors(n_graph, time_ego, ego_lat, ddt, ego_lat_sigma);
    TGraphErrors *graph_lat_ipa = new TGraphErrors(n_gps_ipa, time_ipa, ipa_lat, ddt, ipa_lat_sigma);

    n_graph = 21;
    TH1D *h_dlat_ego = new TH1D("dlat_ego", "dlat (ego)", n_graph, -1.05 * rms_ego, 1.05 * rms_ego);
    TH1D *h_dlat_ipa = new TH1D("dlat_ipa", "dlat (ipa)", n_graph, -1.05 * rms_ipa, 1.05 * rms_ipa);
    TH1D *h_dlat_pulls_ego = new TH1D("dlat_ego_pulls", "dlat pulls (ego)", n_graph, -1.05, 1.05);
    TH1D *h_dlat_pulls_ipa = new TH1D("dlat_ipa_pulls", "dlat pulls (ipa)", n_graph, -1.05, 1.05);
    TH1D *h_dlat_pulls_a_ego = new TH1D("dlat_ego_pulls_a", "dlat pulls (ego)", n_graph, -1.05, 1.05);
    TH1D *h_dlat_pulls_a_ipa = new TH1D("dlat_ipa_pulls_a", "dlat pulls (ipa)", n_graph, -1.05, 1.05);
    for (i_gps = 1; i_gps < n_gps_ego; ++i_gps) {
        double dlat = ego_lat[i_gps] - ego_lat[i_gps - 1];
        double sigma = (ego_lat_sigma[i_gps] + ego_lat_sigma[i_gps-1]) * 0.5;
        h_dlat_ego->Fill(dlat);
        h_dlat_pulls_ego->Fill(dlat / sigma);
        h_dlat_pulls_a_ego->Fill(dlat / rms_ego);
    }

    for (i_gps = 1; i_gps < n_gps_ipa; ++i_gps) {
        double dlat = ipa_lat[i_gps] - ipa_lat[i_gps - 1];
        double sigma = (ipa_lat_sigma[i_gps] + ipa_lat_sigma[i_gps-1]) * 0.5;
        h_dlat_ipa->Fill(dlat);
        h_dlat_pulls_ipa->Fill(dlat / sigma);
        h_dlat_pulls_a_ipa->Fill(dlat / rms_ipa);
    }

    std::string ofile = data_directory + "/analysis.root";
    TFile fp(ofile.c_str(), "recreate");
    graph_dlat->Write("dlat");
    graph_lat_ego->Write("lat_ego");
    graph_lat_ipa->Write("lat_ipa");
    fp.Write();
    fp.Close();

    TCanvas canvas("main", "main", 1200, 800);
    TGraphErrors *graph;

    purdifyPlots(h_dlat_ipa, h_dlat_ego, "#Delta(LAT)", "distance [m]");

    h_dlat_ipa->DrawNormalized("hist");
    h_dlat_ego->DrawNormalized("p same");

    legend = new TLegend(0.7, 0.7, 0.8, 0.8);
    legend->AddEntry(h_dlat_ipa, "IPA", "LF");
    legend->AddEntry(h_dlat_ego, "EGO", "LP");
    legend->Draw();

    canvas.Draw();
    canvas.Update();
    canvas.WaitPrimitive();

    delete legend;

    purdifyPlots(h_dlat_pulls_a_ipa, h_dlat_pulls_a_ego, "#Delta(LAT) / RMS", "#Delta(LAT) / RMS");

    h_dlat_pulls_a_ipa->DrawNormalized("hist");
    h_dlat_pulls_a_ego->DrawNormalized("p same");

    legend = new TLegend(0.7, 0.7, 0.8, 0.8);
    legend->AddEntry(h_dlat_pulls_a_ipa, "IPA", "LF");
    legend->AddEntry(h_dlat_pulls_a_ego, "EGO", "LP");
    legend->Draw();

    canvas.Draw();
    canvas.Update();
    canvas.WaitPrimitive();

    delete legend;

    purdifyPlots(h_dlat_pulls_ipa, h_dlat_pulls_ego, "#Delta(LAT) / #sigma", "#Delta(LAT) / #sigma");

    h_dlat_pulls_ipa->DrawNormalized("hist");
    h_dlat_pulls_ego->DrawNormalized("p same");

    legend = new TLegend(0.7, 0.7, 0.8, 0.8);
    legend->AddEntry(h_dlat_pulls_ipa, "IPA", "LF");
    legend->AddEntry(h_dlat_pulls_ego, "EGO", "LP");
    legend->Draw();

    canvas.Draw();
    canvas.Update();
    canvas.WaitPrimitive();

    delete legend;

    graph = graph_lat_ego;
    graph->SetMinimum(min_lat - max_ego_lat_sigma);
    graph->SetMaximum(max_lat + max_ego_lat_sigma);
    graph->SetMarkerStyle(kFullCircle);
    graph->SetMarkerColor(kBlack);
    graph->SetFillColor(kYellow);
    graph->SetFillStyle(3005);
    graph->Draw("AL4");
//    graph = graph_lat_ego;
//    graph->SetMarkerStyle(kFullTriangleUp);
//    graph->Draw("P");
    canvas.Draw();
    canvas.Update();
    canvas.WaitPrimitive();

    getchar();

    return 0;
}

#include <TRandom3.h>

void runExperiment() {
    double mean = 0.0;
    double sigma = 1.0;
    int i, n = 1000000;
    TRandom3 rndm;
    TH1D *hist = new TH1D("hist", "hist", 100, -5.0, 5.0);
    double x0 = rndm.Gaus(mean, sigma);
    for (i=0;i<n;++i) {
        double x1 = rndm.Gaus(mean, sigma);
        double x = x1 - x0;
        hist->Fill(x);
        x0 = x1;
    }
    TCanvas *c = new TCanvas("test", "test", 800, 600);
    hist->Draw("hist");
    c->Draw();
    c->Update();
    c->WaitPrimitive();
}


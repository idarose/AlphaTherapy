#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TKey.h>
#include <TClass.h>
#include <TIterator.h>
#include <vector>
#include "TROOT.h"
#include <iostream>

// g++ `root-config --cflags` combineHistogramCases.cpp `root-config --libs` -o combine


void CombineHistogramsPC3PIP(std::string cellGeometry, int activity)
{
    // Open the ROOT file
    std::string filename1 = "Output_" + cellGeometry + "/Output_PC3_PIP_" + std::to_string(activity) + "kBq_case_1.root";
    std::string filename2 = "Output_" + cellGeometry + "/Output_PC3_PIP_" + std::to_string(activity) + "kBq_case_2.root";
    std::string filename3 = "Output_" + cellGeometry + "/Output_PC3_PIP_" + std::to_string(activity) + "kBq_case_3.root";
    std::string filename4 = "Output_" + cellGeometry + "/Output_PC3_PIP_" + std::to_string(activity) + "kBq_case_4.root";


    TFile *file1 = TFile::Open(filename1.c_str(), "READ");

    std::vector<TH1D*> vec1DHists_Case1;
    std::vector<TH2D*> vec2DHists_Case1;


    if (file1 && !file1->IsZombie()) {
        TList* list = file1->GetListOfKeys();
        TIter next(list);
        TKey* key;

        // Loop over all keys and retrieve objects
        while ((key = (TKey*)next())) {
            TObject* obj = key->ReadObj();  // Declare 'obj' here
            TClass* cl = gROOT->GetClass(obj->ClassName());
            if (!obj || !cl) {
                delete obj;  // Guard against memory leak if obj is not null
                continue;
            }

            if (cl->InheritsFrom("TH1D")) {
                TH1D* histo = (TH1D*)obj->Clone(); // Clone the histogram
                histo->SetDirectory(0);
                vec1DHists_Case1.push_back(histo);
            }
            else if (cl->InheritsFrom("TH2D")) {
                TH2D* histo2 = (TH2D*)obj->Clone(); // Clone the 2D histogram
                histo2->SetDirectory(0);
                vec2DHists_Case1.push_back(histo2);
            }

            // Ensure you delete the original object only when you are sure you will not use it anymore
            delete obj;
        }

    }
    else {
        std::cerr << "Could not open file or file is corrupted: " << filename1 << std::endl;
    }

    file1->Close();


    TFile *file2 = TFile::Open(filename2.c_str(), "READ");

    std::vector<TH1D*> vec1DHists_Case2;
    std::vector<TH2D*> vec2DHists_Case2;


    if (file2 && !file2->IsZombie()) {
        TList* list = file2->GetListOfKeys();
        TIter next(list);
        TKey* key;

        // Loop over all keys and retrieve objects
        while ((key = (TKey*)next())) {
            TObject* obj = key->ReadObj();  // Declare 'obj' here
            TClass* cl = gROOT->GetClass(obj->ClassName());
            if (!obj || !cl) {
                delete obj;  // Guard against memory leak if obj is not null
                continue;
            }

            if (cl->InheritsFrom("TH1D")) {
                TH1D* histo = (TH1D*)obj->Clone(); // Clone the histogram
                histo->SetDirectory(0);
                vec1DHists_Case2.push_back(histo);
            }
            else if (cl->InheritsFrom("TH2D")) {
                TH2D* histo2 = (TH2D*)obj->Clone(); // Clone the 2D histogram
                histo2->SetDirectory(0);
                vec2DHists_Case2.push_back(histo2);
            }

            // Ensure you delete the original object only when you are sure you will not use it anymore
            delete obj;
        }

    }
    else {
        std::cerr << "Could not open file or file is corrupted: " << filename1 << std::endl;
    }

    file2->Close();

    for(int i=0; i<vec1DHists_Case2.size(); i++)
    {
        (vec1DHists_Case1[i])->Add(vec1DHists_Case2[i]);
    }
    for(int i=0; i<vec2DHists_Case2.size(); i++)
    {
        (vec2DHists_Case1[i])->Add(vec2DHists_Case2[i]);
    }

    vec1DHists_Case2.clear();
    vec2DHists_Case2.clear();


    TFile *file3= TFile::Open(filename3.c_str(), "READ");

    std::vector<TH1D*> vec1DHists_Case3;
    std::vector<TH2D*> vec2DHists_Case3;


    if (file3 && !file3->IsZombie()) {
        TList* list = file3->GetListOfKeys();
        TIter next(list);
        TKey* key;

        // Loop over all keys and retrieve objects
        while ((key = (TKey*)next())) {
            TObject* obj = key->ReadObj();  // Declare 'obj' here
            TClass* cl = gROOT->GetClass(obj->ClassName());
            if (!obj || !cl) {
                delete obj;  // Guard against memory leak if obj is not null
                continue;
            }

            if (cl->InheritsFrom("TH1D")) {
                TH1D* histo = (TH1D*)obj->Clone(); // Clone the histogram
                histo->SetDirectory(0);
                vec1DHists_Case3.push_back(histo);
            }
            else if (cl->InheritsFrom("TH2D")) {
                TH2D* histo2 = (TH2D*)obj->Clone(); // Clone the 2D histogram
                histo2->SetDirectory(0);
                vec2DHists_Case3.push_back(histo2);
            }

            // Ensure you delete the original object only when you are sure you will not use it anymore
            delete obj;
        }

    }
    else {
        std::cerr << "Could not open file or file is corrupted: " << filename1 << std::endl;
    }

    file3->Close();


    for(int i=0; i<vec1DHists_Case3.size(); i++)
    {
        (vec1DHists_Case1[i])->Add(vec1DHists_Case3[i]);
    }
    for(int i=0; i<vec2DHists_Case3.size(); i++)
    {
        (vec2DHists_Case1[i])->Add(vec2DHists_Case3[i]);
    }

    vec1DHists_Case3.clear();
    vec2DHists_Case3.clear();


    TFile *file4= TFile::Open(filename4.c_str(), "READ");

    std::vector<TH1D*> vec1DHists_Case4;
    std::vector<TH2D*> vec2DHists_Case4;


    if (file4 && !file4->IsZombie()) {
        TList* list = file4->GetListOfKeys();
        TIter next(list);
        TKey* key;

        // Loop over all keys and retrieve objects
        while ((key = (TKey*)next())) {
            TObject* obj = key->ReadObj();  // Declare 'obj' here
            TClass* cl = gROOT->GetClass(obj->ClassName());
            if (!obj || !cl) {
                delete obj;  // Guard against memory leak if obj is not null
                continue;
            }

            if (cl->InheritsFrom("TH1D")) {
                TH1D* histo = (TH1D*)obj->Clone(); // Clone the histogram
                histo->SetDirectory(0);
                vec1DHists_Case4.push_back(histo);
            }
            else if (cl->InheritsFrom("TH2D")) {
                TH2D* histo2 = (TH2D*)obj->Clone(); // Clone the 2D histogram
                histo2->SetDirectory(0);
                vec2DHists_Case4.push_back(histo2);
            }

            // Ensure you delete the original object only when you are sure you will not use it anymore
            delete obj;
        }

    }
    else {
        std::cerr << "Could not open file or file is corrupted: " << filename1 << std::endl;
    }

    file4->Close();


    for(int i=0; i<vec1DHists_Case4.size(); i++)
    {
        (vec1DHists_Case1[i])->Add(vec1DHists_Case3[i]);
    }
    for(int i=0; i<vec2DHists_Case4.size(); i++)
    {
        (vec2DHists_Case1[i])->Add(vec2DHists_Case4[i]);
    }

    vec1DHists_Case4.clear();
    vec2DHists_Case4.clear();



    std::string outputName = "Output_" + cellGeometry + "/Output_PC3_PIP_" + std::to_string(activity) + "kBq.root";
    auto output = new TFile(outputName.c_str(), "RECREATE");

    for(auto & entry : vec1DHists_Case1)
    {
        entry->Scale(1./4.);
        entry->Write();
    }
    output->Write();
    output->Close();

}

void CombineHistogramsPC3Flu(std::string cellGeometry, int activity)
{
    // Open the ROOT file
    std::string filename1 = "Output_" + cellGeometry + "/Output_PC3_Flu_" + std::to_string(activity) + "kBq_case_1.root";
    std::string filename2 = "Output_" + cellGeometry + "/Output_PC3_Flu_" + std::to_string(activity) + "kBq_case_2.root";


    TFile *file1 = TFile::Open(filename1.c_str(), "READ");

    std::vector<TH1D*> vec1DHists_Case1;
    std::vector<TH2D*> vec2DHists_Case1;


    if (file1 && !file1->IsZombie()) {
        TList* list = file1->GetListOfKeys();
        TIter next(list);
        TKey* key;

        // Loop over all keys and retrieve objects
        while ((key = (TKey*)next())) {
            TObject* obj = key->ReadObj();  // Declare 'obj' here
            TClass* cl = gROOT->GetClass(obj->ClassName());
            if (!obj || !cl) {
                delete obj;  // Guard against memory leak if obj is not null
                continue;
            }

            if (cl->InheritsFrom("TH1D")) {
                TH1D* histo = (TH1D*)obj->Clone(); // Clone the histogram
                histo->SetDirectory(0);
                vec1DHists_Case1.push_back(histo);
            }
            else if (cl->InheritsFrom("TH2D")) {
                TH2D* histo2 = (TH2D*)obj->Clone(); // Clone the 2D histogram
                histo2->SetDirectory(0);
                vec2DHists_Case1.push_back(histo2);
            }

            // Ensure you delete the original object only when you are sure you will not use it anymore
            delete obj;
        }

    }
    else {
        std::cerr << "Could not open file or file is corrupted: " << filename1 << std::endl;
    }

    file1->Close();


    TFile *file2 = TFile::Open(filename2.c_str(), "READ");

    std::vector<TH1D*> vec1DHists_Case2;
    std::vector<TH2D*> vec2DHists_Case2;


    if (file2 && !file2->IsZombie()) {
        TList* list = file2->GetListOfKeys();
        TIter next(list);
        TKey* key;

        // Loop over all keys and retrieve objects
        while ((key = (TKey*)next())) {
            TObject* obj = key->ReadObj();  // Declare 'obj' here
            TClass* cl = gROOT->GetClass(obj->ClassName());
            if (!obj || !cl) {
                delete obj;  // Guard against memory leak if obj is not null
                continue;
            }

            if (cl->InheritsFrom("TH1D")) {
                TH1D* histo = (TH1D*)obj->Clone(); // Clone the histogram
                histo->SetDirectory(0);
                vec1DHists_Case2.push_back(histo);
            }
            else if (cl->InheritsFrom("TH2D")) {
                TH2D* histo2 = (TH2D*)obj->Clone(); // Clone the 2D histogram
                histo2->SetDirectory(0);
                vec2DHists_Case2.push_back(histo2);
            }

            // Ensure you delete the original object only when you are sure you will not use it anymore
            delete obj;
        }

    }
    else {
        std::cerr << "Could not open file or file is corrupted: " << filename1 << std::endl;
    }

    file2->Close();

    for(int i=0; i<vec1DHists_Case2.size(); i++)
    {
        (vec1DHists_Case1[i])->Add(vec1DHists_Case2[i]);
    }
    for(int i=0; i<vec2DHists_Case2.size(); i++)
    {
        (vec2DHists_Case1[i])->Add(vec2DHists_Case2[i]);
    }

    vec1DHists_Case2.clear();
    vec2DHists_Case2.clear();


    std::string outputName = "Output_" + cellGeometry + "/Output_PC3_Flu_" + std::to_string(activity) + "kBq.root";
    auto output = new TFile(outputName.c_str(), "RECREATE");

    for(auto & entry : vec1DHists_Case1)
    {
        entry->Scale(1./2.);
        entry->Write();
    }
    output->Write();
    output->Close();

}

void CombineHistogramsC42(std::string cellGeometry, int activity)
{
    // Open the ROOT file
    std::string filename1 = "Output_" + cellGeometry + "/Output_C4_2_" + std::to_string(activity) + "kBq_case_1.root";


    TFile *file1 = TFile::Open(filename1.c_str(), "READ");

    std::vector<TH1D*> vec1DHists_Case1;


    if (file1 && !file1->IsZombie()) {
        TList* list = file1->GetListOfKeys();
        TIter next(list);
        TKey* key;

        // Loop over all keys and retrieve objects
        while ((key = (TKey*)next())) {
            TObject* obj = key->ReadObj();  // Declare 'obj' here
            TClass* cl = gROOT->GetClass(obj->ClassName());
            if (!obj || !cl) {
                delete obj;  // Guard against memory leak if obj is not null
                continue;
            }

            if (cl->InheritsFrom("TH1D")) {
                TH1D* histo = (TH1D*)obj->Clone(); // Clone the histogram
                histo->SetDirectory(0);
                vec1DHists_Case1.push_back(histo);
            }
            else if (cl->InheritsFrom("TH2D")) {
                TH2D* histo2 = (TH2D*)obj->Clone(); // Clone the 2D histogram
                histo2->SetDirectory(0);
                vec2DHists_Case1.push_back(histo2);
            }

            // Ensure you delete the original object only when you are sure you will not use it anymore
            delete obj;
        }

    }
    else {
        std::cerr << "Could not open file or file is corrupted: " << filename1 << std::endl;
    }

    file1->Close();


    std::string outputName = "Output_" + cellGeometry + "/Output_C4_2_" + std::to_string(activity) + "kBq.root";
    auto output = new TFile(outputName.c_str(), "RECREATE");

    for(auto & entry : vec1DHists_Case1)
    {
        entry->Scale(1./1.);
        entry->Write();
    }
    output->Write();
    output->Close();
}

int main(int argc, char *argv[]) {

    //------------------–----------
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " string_argument integer_argument" << std::endl;
        return 1;
    }

    //------------------–----------
    std::string cellGeometry = argv[1];
    int activity;


    std::vector<int> validActivities = {1,3,5,10,25,50,75,100,150};

    //--------------------------
    try {
        activity = std::stoi(argv[2]);
        bool activityIsValid = false;
        for(auto & entry : validActivities){
            if(entry == activity)
            {
                activityIsValid = true;
            }
        }
        if(!activityIsValid)
        {
            throw std::invalid_argument("Activity case not valid.");
        }
    } catch (const std::invalid_argument& e) {
        std::cerr << "Error: second argument is not a valid integer! " << e.what() << std::endl;
        return 2;
    } catch (const std::out_of_range& e) {
        std::cerr << "Error: second argument is out of range for an integer! " << e.what() << std::endl;
        return 3;
    }

    // Call the function with the name of your ROOT file
    CombineHistogramsPC3PIP(cellGeometry, activity);
    CombineHistogramsPC3Flu(cellGeometry, activity);
    CombineHistogramsC42(cellGeometry, activity);

    return 0;
}
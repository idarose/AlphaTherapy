#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TKey.h>
#include <TClass.h>
#include <TIterator.h>
#include <vector>

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

int main() {
    // Call the function with the name of your ROOT file
    CombineHistogramsPC3PIP("D12RP", 1);
    CombineHistogramsPC3PIP("D12RP", 3);
    CombineHistogramsPC3PIP("D12RP", 5);
    CombineHistogramsPC3PIP("D12RP", 10);
    CombineHistogramsPC3PIP("D12RP", 25);
    CombineHistogramsPC3PIP("D12RP", 50);
    CombineHistogramsPC3PIP("D12RP", 75);
    CombineHistogramsPC3PIP("D12RP", 100);
    CombineHistogramsPC3PIP("D12RP", 150);
    return 0;
}
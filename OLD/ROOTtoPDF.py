void print_all_canvases(const char* filename = "myfile.root", const char* outpdf = "output.pdf") {
  // Open the ROOT file
  TFile* f = TFile::Open(filename);
  if (!f || f->IsZombie()) {
    Error("print_all_canvases", "Cannot open file: %s", filename);
    return;
  }

  // Get list of keys in the ROOT file
  TIter next(f->GetListOfKeys());
  TKey* key;

  // Start multipage PDF
  TString pdfname = outpdf;
  TString pdfopen = pdfname + "[";
  TString pdfclose = pdfname + "]";

  // Dummy canvas to start/close PDF
  TCanvas* dummy = new TCanvas();
  dummy->Print(pdfopen); // Start PDF

  // Loop over keys and print all canvases
  while ((key = (TKey*)next())) {
    TObject* obj = key->ReadObj();

    if (obj->InheritsFrom("TCanvas")) {
      TCanvas* c = (TCanvas*)obj;
      c->Draw();
      c->Print(pdfname);
    }
  }

  dummy->Print(pdfclose); // Close PDF
  delete dummy;

  f->Close();
  delete f;

  printf("✅ Canvases saved to: %s\n", outpdf);
}

// To use:
//   root
//   .L ROOTtoPDF.C
//   ROOTtoPDF("moveme.root");
//   ROOTtoPDF("moveme.root","custom.pdf");

void ROOTtoPDF(const char* filename = "moveme.root",
               const char* outpdf   = "")
{
  // ----------------------------------
  // Build output PDF name if not given
  // ----------------------------------
  TString pdfname;
  if (strlen(outpdf) == 0) {
    pdfname = filename;
    if (!pdfname.EndsWith(".root"))
      pdfname += ".root";
    pdfname.ReplaceAll(".root", ".pdf");
  } else {
    pdfname = outpdf;
  }

  // ----------------------------------
  // Open the ROOT file
  // ----------------------------------
  TFile* f = TFile::Open(filename);
  if (!f || f->IsZombie()) {
    Error("ROOTtoPDF", "Cannot open file: %s", filename);
    return;
  }

  // ----------------------------------
  // Prepare multipage PDF
  // ----------------------------------
  TString pdfopen  = pdfname + "[";
  TString pdfclose = pdfname + "]";

  TCanvas* dummy = new TCanvas();
  dummy->Print(pdfopen);

  // ----------------------------------
  // Loop over keys and print canvases
  // ----------------------------------
  TIter next(f->GetListOfKeys());
  TKey* key;

  while ((key = (TKey*)next())) {
    TObject* obj = key->ReadObj();
    if (obj->InheritsFrom("TCanvas")) {
      TCanvas* c = (TCanvas*)obj;
      c->Draw();
      c->Print(pdfname);
    }
  }

  dummy->Print(pdfclose);

  delete dummy;
  f->Close();
  delete f;

  printf("✅ Canvases saved to: %s\n", pdfname.Data());
}

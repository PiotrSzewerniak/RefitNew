void dohist()
{
	gROOT->ProcessLine(".L hist.C");
	gROOT->ProcessLine("hist t");
	gROOT->ProcessLine("t.Loop()");
}
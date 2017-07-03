﻿using NUnit.Framework;
using ProteoformSuiteInternal;
using Proteomics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace Test
{
    [TestFixture]
    class TestGoTermNumber
    {
        [Test]
        public void testLogOddsRatio()
        {
            int q = 1;//number of enriched proteins with the term
            int k = 2;//number of enriched proteins
            int m = 2;//number of proteins in the background with the term
            int t = 4;//number of proteins in the background

            DatabaseReference d = new DatabaseReference("GO", ":1", new List<Tuple<string, string>> { new Tuple<string, string>("term", "P:one") });
            GoTerm g = new GoTerm(d);
            GoTermNumber gtn = new GoTermNumber(g, q, k, m, t);
            Assert.AreEqual(0, gtn.log_odds_ratio);
            //Assert.AreEqual(0.833333m, Math.Round((decimal)gtn.p_value.Truncate(7), 6));
            Assert.AreEqual(0.833333m, Math.Round(gtn.p_value, 6));

            q = 1;//number of enriched proteins with the term
            k = 2;//number of enriched proteins
            m = 4;//number of proteins in the background with the term
            t = 4;//number of proteins in the background

            gtn = new GoTermNumber(g, q, k, m, t);
            Assert.AreEqual(-1, gtn.log_odds_ratio);
            //Assert.AreEqual(1m, Math.Round((decimal)gtn.p_value.Truncate(7), 6));
            Assert.AreEqual(1m, Math.Round(gtn.p_value, 6));

            q = 2;//number of enriched proteins with the term
            k = 2;//number of enriched proteins
            m = 2;//number of proteins in the background with the term
            t = 4;//number of proteins in the background

            gtn = new GoTermNumber(g, q, k, m, t);
            Assert.AreEqual(1, gtn.log_odds_ratio);
            //Assert.AreEqual(0.166667m, Math.Round((decimal)gtn.p_value.Truncate(7), 6));
            Assert.AreEqual(0.166667m, Math.Round(gtn.p_value, 6));
        }

        [Test]
        public void testBinomialCoefficient()
        {
            int n = 10;
            int k = 2;
            Assert.AreEqual(45, (int)GoTermNumber.binomialCoefficient(n, k));
        }

        [Test]
        public void testGoTerm_pValue()
        {
            int q = 2; //count of proteins in selected subset with the particular Go term
            int k = 4; //count of proteins in selected subset
            int m = 4; //count of proteins in background with the particular Go term
            int t = 10; //count of proteins in background

            Assert.AreEqual(0.54762, Math.Round((double)GoTermNumber.GoTerm_pValue(q, k, m, t), 5));
        }

        [Test]
        public void get_interesting_goterm_families()
        {
            SaveState.lollipop = new Lollipop();
            DatabaseReference d1 = new DatabaseReference("GO", "GO:1", new List<Tuple<string, string>> { new Tuple<string, string>("term", "P:1") });
            DatabaseReference d2 = new DatabaseReference("GO", "GO:2", new List<Tuple<string, string>> { new Tuple<string, string>("term", "P:2") });
            DatabaseReference d3 = new DatabaseReference("GO", "GO:1", new List<Tuple<string, string>> { new Tuple<string, string>("term", "P:1") });
            GoTerm g1 = new GoTerm(d1);
            GoTerm g2 = new GoTerm(d2);
            GoTerm g3 = new GoTerm(d3);
            ProteinWithGoTerms p1 = new ProteinWithGoTerms("", "T1", new List<Tuple<string, string>> { new Tuple<string, string>("", "") }, new Dictionary<int, List<Modification>>(), new int?[] { 0 }, new int?[] { 0 }, new string[] { "" }, "T2", "T3", true, false, new List<DatabaseReference> { d1 }, new List<GoTerm> { g1 }, new List<DisulfideBond>());
            ProteinWithGoTerms p2 = new ProteinWithGoTerms("", "T2", new List<Tuple<string, string>> { new Tuple<string, string>("", "") }, new Dictionary<int, List<Modification>>(), new int?[] { 0 }, new int?[] { 0 }, new string[] { "" }, "T2", "T3", true, false, new List<DatabaseReference> { d2 }, new List<GoTerm> { g2 }, new List<DisulfideBond>());
            ProteinWithGoTerms p3 = new ProteinWithGoTerms("", "T3", new List<Tuple<string, string>> { new Tuple<string, string>("", "") }, new Dictionary<int, List<Modification>>(), new int?[] { 0 }, new int?[] { 0 }, new string[] { "" }, "T2", "T3", true, false, new List<DatabaseReference> { d3 }, new List<GoTerm> { g3 }, new List<DisulfideBond>());
            Dictionary<InputFile, Protein[]> dict = new Dictionary<InputFile, Protein[]>
            {
                { new InputFile("fake.txt", Purpose.ProteinDatabase), new Protein[] { p1 } },
                { new InputFile("fake.txt", Purpose.ProteinDatabase), new Protein[] { p2 } },
                { new InputFile("fake.txt", Purpose.ProteinDatabase), new Protein[] { p3 } },
            };
            ExperimentalProteoform e1 = ConstructorsForTesting.ExperimentalProteoform("E");
            ExperimentalProteoform e2 = ConstructorsForTesting.ExperimentalProteoform("E");
            e1.quant.intensitySum = 1;
            e1.quant.FDR = 0;
            e1.quant.logFoldChange = 1;
            e2.quant.intensitySum = 1;
            e2.quant.FDR = 0;
            e2.quant.logFoldChange = 1;
            TheoreticalProteoform t = ConstructorsForTesting.make_a_theoretical("T1_T1_asdf", p1, dict);
            TheoreticalProteoform u = ConstructorsForTesting.make_a_theoretical("T2_T1_asdf_asdf", p2, dict);
            TheoreticalProteoform v = ConstructorsForTesting.make_a_theoretical("T3_T1_asdf_Asdf_Asdf", p3, dict);
            t.ExpandedProteinList = new List<ProteinWithGoTerms> { p1 };
            u.ExpandedProteinList = new List<ProteinWithGoTerms> { p2 };
            v.ExpandedProteinList = new List<ProteinWithGoTerms> { p3 };
            make_relation(e1, t);
            //make_relation(e1, v); // we don't allow this to happen anymore... we only allow one ET conntection per E
            make_relation(e2, u);
            ProteoformFamily f = new ProteoformFamily(e1); // two theoreticals with the same GoTerms... expecting one GoTerm number but two theoretical proteins (now only one)
            ProteoformFamily h = new ProteoformFamily(e2);
            f.construct_family();
            f.identify_experimentals();
            h.construct_family();
            h.identify_experimentals();
            List<ProteoformFamily> families = new List<ProteoformFamily> { f, h };
            t.family = f;
            v.family = f;
            e1.family = f;
            u.family = h;
            e2.family = h;
            List<ExperimentalProteoform> fake_significant = new List<ExperimentalProteoform> { e1 };
            List<ProteinWithGoTerms> significant_proteins = SaveState.lollipop.getInducedOrRepressedProteins(fake_significant, 0, 1, 0);
            List<GoTermNumber> gtn = SaveState.lollipop.getGoTermNumbers(significant_proteins, new List<ProteinWithGoTerms> { p1, p2, p3 });
            Assert.AreEqual(1, significant_proteins.Count);
            Assert.AreEqual(1, gtn.Count);
            Assert.AreEqual("1", gtn.First().Id);
            Assert.AreEqual(0 - (decimal)Math.Log(2d / 3d, 2), gtn.First().log_odds_ratio);

            List<ProteoformFamily> fams = SaveState.lollipop.getInterestingFamilies(gtn, families);
            Assert.AreEqual(1, fams.Count);
            Assert.AreEqual(1, fams[0].theoretical_proteoforms.Count);
        }

        [Test]
        public void test_goterm_analysis()
        {
            SaveState.lollipop = new Lollipop();
            DatabaseReference d1 = new DatabaseReference("GO", "GO:1", new List<Tuple<string, string>> { new Tuple<string, string>("term", "P:1") });
            DatabaseReference d2 = new DatabaseReference("GO", "GO:2", new List<Tuple<string, string>> { new Tuple<string, string>("term", "P:2") });
            DatabaseReference d3 = new DatabaseReference("GO", "GO:1", new List<Tuple<string, string>> { new Tuple<string, string>("term", "P:1") });
            GoTerm g1 = new GoTerm(d1);
            GoTerm g2 = new GoTerm(d2);
            GoTerm g3 = new GoTerm(d3);
            ProteinWithGoTerms p1 = new ProteinWithGoTerms("", "T1", new List<Tuple<string, string>> { new Tuple<string, string>("", "") }, new Dictionary<int, List<Modification>>(), new int?[] { 0 }, new int?[] { 0 }, new string[] { "" }, "T2", "T3", true, false, new List<DatabaseReference> { d1 }, new List<GoTerm> { g1 }, new List<DisulfideBond>());
            ProteinWithGoTerms p2 = new ProteinWithGoTerms("", "T2", new List<Tuple<string, string>> { new Tuple<string, string>("", "") }, new Dictionary<int, List<Modification>>(), new int?[] { 0 }, new int?[] { 0 }, new string[] { "" }, "T2", "T3", true, false, new List<DatabaseReference> { d2 }, new List<GoTerm> { g2 }, new List<DisulfideBond>());
            ProteinWithGoTerms p3 = new ProteinWithGoTerms("", "T3", new List<Tuple<string, string>> { new Tuple<string, string>("", "") }, new Dictionary<int, List<Modification>>(), new int?[] { 0 }, new int?[] { 0 }, new string[] { "" }, "T2", "T3", true, false, new List<DatabaseReference> { d3 }, new List<GoTerm> { g3 }, new List<DisulfideBond>());
            Dictionary<InputFile, Protein[]> dict = new Dictionary<InputFile, Protein[]> {
                { new InputFile("fake.txt", Purpose.ProteinDatabase), new Protein[] { p1 } },
                { new InputFile("fake.txt", Purpose.ProteinDatabase), new Protein[] { p2 } },
                { new InputFile("fake.txt", Purpose.ProteinDatabase), new Protein[] { p3 } },
            };
            ExperimentalProteoform e1 = ConstructorsForTesting.ExperimentalProteoform("E");
            ExperimentalProteoform e2 = ConstructorsForTesting.ExperimentalProteoform("E");
            e1.quant.intensitySum = 1;
            e1.quant.FDR = 0;
            e1.quant.logFoldChange = 1;
            e2.quant.intensitySum = 1;
            e2.quant.FDR = 0;
            e2.quant.logFoldChange = 1;
            TheoreticalProteoform t = ConstructorsForTesting.make_a_theoretical("T1_T1_asdf", p1, dict);
            TheoreticalProteoform u = ConstructorsForTesting.make_a_theoretical("T2_T1_asdf_asdf", p2, dict);
            TheoreticalProteoform v = ConstructorsForTesting.make_a_theoretical("T3_T1_asdf_Asdf_Asdf", p3, dict);
            t.ExpandedProteinList = new List<ProteinWithGoTerms> { p1 };
            u.ExpandedProteinList = new List<ProteinWithGoTerms> { p2 };
            v.ExpandedProteinList = new List<ProteinWithGoTerms> { p3 };
            make_relation(e1, t);
            make_relation(e1, v);
            make_relation(e2, u);
            ProteoformFamily f = new ProteoformFamily(e1); // two theoreticals with the same GoTerms... expecting one GoTerm number but two theoretical proteins
            ProteoformFamily h = new ProteoformFamily(e2);
            f.construct_family();
            f.identify_experimentals();
            h.construct_family();
            h.identify_experimentals();
            List<ProteoformFamily> families = new List<ProteoformFamily> { f, h };
            t.family = f;
            v.family = f;
            e1.family = f;
            u.family = h;
            e2.family = h;
            SaveState.lollipop.inducedOrRepressedProteins = SaveState.lollipop.getInducedOrRepressedProteins(new List<ExperimentalProteoform> { e1 }, 0, 1, 0);
            SaveState.lollipop.allTheoreticalProteins = true;
            SaveState.lollipop.theoretical_database.expanded_proteins = new ProteinWithGoTerms[] { p1, p2, p3 };
            SaveState.lollipop.GO_analysis();
            Assert.AreEqual(1, SaveState.lollipop.inducedOrRepressedProteins.Count);  // only taking one ET connection by definition in forming ET relations; only one is used in identify theoreticals
            Assert.AreEqual(1, SaveState.lollipop.goTermNumbers.Count);
            Assert.AreEqual("1", SaveState.lollipop.goTermNumbers.First().Id);
            Assert.AreEqual(0 - (decimal)Math.Log(2d / 3d, 2), SaveState.lollipop.goTermNumbers.First().log_odds_ratio);

            List<ProteoformFamily> fams = SaveState.lollipop.getInterestingFamilies(SaveState.lollipop.goTermNumbers, families);
            Assert.AreEqual(1, fams.Count);
            Assert.AreEqual(2, fams[0].theoretical_proteoforms.Count);
        }

        public void make_relation(Proteoform p1, Proteoform p2)
        {
            ProteoformRelation pp = new ProteoformRelation(p1, p2, ProteoformComparison.ExperimentalExperimental, 0, TestContext.CurrentContext.TestDirectory);
            DeltaMassPeak ppp = new DeltaMassPeak(pp, new List<ProteoformRelation> { pp });
            pp.Accepted = true;
            pp.peak = ppp;
            pp.Accepted = true;
            ppp.Accepted = true;
            p1.relationships.Add(pp);
            p2.relationships.Add(pp);
        }

        [Test]
        public void test_goterm_analysis_with_custom_list()
        {
            SaveState.lollipop = new Lollipop();
            DatabaseReference d1 = new DatabaseReference("GO", "GO:1", new List<Tuple<string, string>> { new Tuple<string, string>("term", "P:1") });
            DatabaseReference d2 = new DatabaseReference("GO", "GO:2", new List<Tuple<string, string>> { new Tuple<string, string>("term", "P:2") });
            DatabaseReference d3 = new DatabaseReference("GO", "GO:1", new List<Tuple<string, string>> { new Tuple<string, string>("term", "P:1") });
            GoTerm g1 = new GoTerm(d1);
            GoTerm g2 = new GoTerm(d2);
            GoTerm g3 = new GoTerm(d3);
            ProteinWithGoTerms p1 = new ProteinWithGoTerms("", "T1", new List<Tuple<string, string>> { new Tuple<string, string>("", "") }, new Dictionary<int, List<Modification>>(), new int?[] { 0 }, new int?[] { 0 }, new string[] { "" }, "T2", "T3", true, false, new List<DatabaseReference> { d1 }, new List<GoTerm> { g1 }, new List<DisulfideBond>());
            ProteinWithGoTerms p2 = new ProteinWithGoTerms("", "T2", new List<Tuple<string, string>> { new Tuple<string, string>("", "") }, new Dictionary<int, List<Modification>>(), new int?[] { 0 }, new int?[] { 0 }, new string[] { "" }, "T2", "T3", true, false, new List<DatabaseReference> { d2 }, new List<GoTerm> { g2 }, new List<DisulfideBond>());
            ProteinWithGoTerms p3 = new ProteinWithGoTerms("", "T3", new List<Tuple<string, string>> { new Tuple<string, string>("", "") }, new Dictionary<int, List<Modification>>(), new int?[] { 0 }, new int?[] { 0 }, new string[] { "" }, "T2", "T3", true, false, new List<DatabaseReference> { d3 }, new List<GoTerm> { g3 }, new List<DisulfideBond>());
            Dictionary<InputFile, Protein[]> dict = new Dictionary<InputFile, Protein[]> {
                { new InputFile("fake.txt", Purpose.ProteinDatabase), new Protein[] { p1 } },
                { new InputFile("fake.txt", Purpose.ProteinDatabase), new Protein[] { p2 } },
                { new InputFile("fake.txt", Purpose.ProteinDatabase), new Protein[] { p3 } },
            };
            ExperimentalProteoform e1 = ConstructorsForTesting.ExperimentalProteoform("E");
            ExperimentalProteoform e2 = ConstructorsForTesting.ExperimentalProteoform("E");
            e1.quant.intensitySum = 1;
            e1.quant.FDR = 0;
            e1.quant.logFoldChange = 1;
            e2.quant.intensitySum = 1;
            e2.quant.FDR = 0;
            e2.quant.logFoldChange = 1;
            TheoreticalProteoform t = ConstructorsForTesting.make_a_theoretical("T1_T1_asdf", p1, dict);
            TheoreticalProteoform u = ConstructorsForTesting.make_a_theoretical("T2_T1_asdf_asdf", p2, dict);
            TheoreticalProteoform v = ConstructorsForTesting.make_a_theoretical("T3_T1_asdf_Asdf_Asdf", p3, dict);
            t.ExpandedProteinList = new List<ProteinWithGoTerms> { p1 };
            u.ExpandedProteinList = new List<ProteinWithGoTerms> { p2 };
            v.ExpandedProteinList = new List<ProteinWithGoTerms> { p3 };
            make_relation(e1, t);
            make_relation(e1, v);
            make_relation(e2, u);
            ProteoformFamily f = new ProteoformFamily(e1); // two theoreticals with the same GoTerms... expecting one GoTerm number but two theoretical proteins
            ProteoformFamily h = new ProteoformFamily(e2);
            f.construct_family();
            f.identify_experimentals();
            h.construct_family();
            h.identify_experimentals();
            List<ProteoformFamily> families = new List<ProteoformFamily> { f, h };
            t.family = f;
            v.family = f;
            e1.family = f;
            u.family = h;
            e2.family = h;
            SaveState.lollipop.inducedOrRepressedProteins = SaveState.lollipop.getInducedOrRepressedProteins(new List<ExperimentalProteoform> { e1 }, 0, 1, 0);
            SaveState.lollipop.allTheoreticalProteins = true;
            SaveState.lollipop.theoretical_database.expanded_proteins = new ProteinWithGoTerms[] { p1, p2, p3 };
            SaveState.lollipop.backgroundProteinsList = Path.Combine(TestContext.CurrentContext.TestDirectory, "test_protein_list.txt");
            SaveState.lollipop.GO_analysis();
            Assert.AreEqual(1, SaveState.lollipop.inducedOrRepressedProteins.Count);  // only taking one ET connection by definition in forming ET relations; only one is used in identify theoreticals
            Assert.AreEqual(1, SaveState.lollipop.goTermNumbers.Count);
            Assert.AreEqual("1", SaveState.lollipop.goTermNumbers.First().Id);
            Assert.AreEqual(0 - (decimal)Math.Log(2d / 3d, 2), SaveState.lollipop.goTermNumbers.First().log_odds_ratio);

            List<ProteoformFamily> fams = SaveState.lollipop.getInterestingFamilies(SaveState.lollipop.goTermNumbers, families);
            Assert.AreEqual(1, fams.Count);
            Assert.AreEqual(2, fams[0].theoretical_proteoforms.Count);
        }
    }
}

﻿using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;

namespace PS_0._00
{
    public partial class AggregatedProteoforms : Form
    {
        DataTableHandler dataTableHandler = new DataTableHandler();

        public AggregatedProteoforms()
        {
            InitializeComponent();
        }

        private void AggregatedProteoforms_Load(object sender, EventArgs e)
        {
            InitializeSettings();
            dataTableHandler.MassDecimals = 8; //Round to 8 decimal places for mass in this form
            GlobalData.acceptableNeuCodeLightProteoforms = FillAcceptableNeuCodeLightProteoformsDataTable();
            GlobalData.aggregatedProteoforms = CreateAggregatedProteoformsDataTable();
            AggregateNeuCodeLightProteoforms();
        }

        private void FillAggregatesTable()
        {
            string[] rt_column_names = new string[] { "Aggregated Retention Time" };
            string[] intensity_column_names = new string[] { "Aggregated Intensity" };
            string[] mass_column_names = new string[] { "Aggregated Mass" };
            string[] abundance_column_names = new string[] { };
            BindingSource bs_aggregatedProteoforms = dataTableHandler.DisplayWithRoundedDoubles(dgv_AggregatedProteoforms, GlobalData.aggregatedProteoforms, 
                rt_column_names, intensity_column_names, abundance_column_names, mass_column_names);
        }

        private void RoundDoubleColumn(DataTable table, string column_name, int num_decimal_places)
        {
            table.AsEnumerable().ToList().ForEach(p => p.SetField<Double>(column_name, Math.Round(p.Field<Double>(column_name), num_decimal_places)));
        }

        private void InitializeSettings()
        {
            nUP_mass_tolerance.Minimum = 0;
            nUP_mass_tolerance.Maximum = 10;
            nUP_mass_tolerance.Value = 3;

            nUD_RetTimeToleranace.Minimum = 0;
            nUD_RetTimeToleranace.Maximum = 10;
            nUD_RetTimeToleranace.Value = 3;

            nUD_Missed_Monos.Minimum = 0;
            nUD_Missed_Monos.Maximum = 5;
            nUD_Missed_Monos.Value = 3;

            nUD_Missed_Ks.Minimum = 0;
            nUD_Missed_Ks.Maximum = 3;
            nUD_Missed_Ks.Value = 1;

        }

        private DataTable CreateAggregatedProteoformsDataTable()
        {
            DataTable dt = new DataTable();

            dt.Columns.Add("Aggregated Mass", typeof(double));
            dt.Columns.Add("Aggregated Intensity", typeof(double));
            dt.Columns.Add("Aggregated Retention Time", typeof(double));
            dt.Columns.Add("Lysine Count", typeof(int));

            return dt;
        }

        private DataTable FillAcceptableNeuCodeLightProteoformsDataTable()
        {
            //MessageBox.Show("Filling ltNCProteoforms table.");
            DataTable acceptableLtProteoforms = new DataTable();
            acceptableLtProteoforms.Columns.Add("Light Filename", typeof(string));
            acceptableLtProteoforms.Columns.Add("Light No#", typeof(int));
            acceptableLtProteoforms.Columns.Add("Light Mass", typeof(double));
            acceptableLtProteoforms.Columns.Add("Light Mass Corrected", typeof(double));
            acceptableLtProteoforms.Columns.Add("Light Intensity", typeof(double));
            acceptableLtProteoforms.Columns.Add("Light Retention Time", typeof(double));
            acceptableLtProteoforms.Columns.Add("Lysine Count", typeof(int));
            acceptableLtProteoforms.Columns.Add("Aggregated Mass", typeof(double));
            acceptableLtProteoforms.Columns.Add("Aggregated Intensity", typeof(double));
            acceptableLtProteoforms.Columns.Add("Aggregated Retention Time", typeof(double));

            foreach (DataRow row in GlobalData.rawNeuCodePairs.Rows)
            {
                if (bool.Parse(row["Acceptable"].ToString()))
                {
                    string lightFilename = row["Light Filename"].ToString();
                    int lightNumber = int.Parse(row["Light No#"].ToString());
                    double ltMass = double.Parse(row["Light Mass"].ToString());
                    double ltMassCorrected = double.Parse(row["Light Mass Corrected"].ToString());
                    double ltIntensity = double.Parse(row["Light Intensity"].ToString());
                    double ltRetentionTime = double.Parse(row["Apex RT"].ToString());
                    int lysineCount = int.Parse(row["Lysine Count"].ToString());

                    double aggregatedMass = 0;
                    double aggregatedIntensity = 0;
                    double aggreagedRetentionTime = 0;

                    acceptableLtProteoforms.Rows.Add(lightFilename, lightNumber, ltMass, ltMassCorrected, ltIntensity, ltRetentionTime, lysineCount, 
                        aggregatedMass, aggregatedIntensity, aggreagedRetentionTime);
                }
            }
            return acceptableLtProteoforms;
        }

        private void ZeroAggregateMasses()
        {
            foreach (DataRow row in GlobalData.acceptableNeuCodeLightProteoforms.Rows)
            {
                row["Aggregated Mass"] = 0;
            }
        }

        private void AggregateNeuCodeLightProteoforms()
        {
            if (GlobalData.aggregatedProteoforms.Rows.Count > 0) { GlobalData.aggregatedProteoforms.Clear(); }
            while (GlobalData.acceptableNeuCodeLightProteoforms.Select("[Aggregated Mass] = 0").Length > 0)
            {
                DataRow[] zeros = GlobalData.acceptableNeuCodeLightProteoforms.Select("[Aggregated Mass] = 0");
                DataRow maxRow = zeros[0];

                foreach (DataRow row in zeros)
                {
                    bool isMaxIntensity = double.Parse(row["Light Intensity"].ToString()) > double.Parse(maxRow["Light Intensity"].ToString());
                    if (isMaxIntensity) { maxRow = row; }
                }

                double mass = double.Parse(maxRow["Light Mass Corrected"].ToString());
                double retTime = double.Parse(maxRow["Light Retention Time"].ToString());
                int lysineCount = int.Parse(maxRow["Lysine Count"].ToString());

                //The section below allows for missed monoisotopics to be aggregated
                string expression = "(";
                for (int shiftNum = -(Convert.ToInt32(nUD_Missed_Monos.Value)); shiftNum <= (Convert.ToInt32(nUD_Missed_Monos.Value)); shiftNum++)
                {
                    double shift = shiftNum * 1.0015;
                    double low = mass + shift - (mass + shift) / 1000000 * Convert.ToInt32(nUP_mass_tolerance.Value);
                    double high = mass + shift + (mass + shift) / 1000000 * Convert.ToInt32(nUP_mass_tolerance.Value);
                    expression = expression + "[Light Mass Corrected] >= " + low +
                    " and [Light Mass Corrected] <= " + high;
                    if (shiftNum < (Convert.ToInt32(nUD_Missed_Monos.Value)))
                    {
                        expression = expression + " or ";
                    }
                    else
                    {
                        expression = expression + ") and ";
                    }
                }

                //The section below allows for missed lysine counts to be aggregated
                expression = expression + "(";
                for (int shiftNum = -(Convert.ToInt32(nUD_Missed_Ks.Value)); shiftNum <= (Convert.ToInt32(nUD_Missed_Ks.Value)); shiftNum++)
                {
                    int K_to_Add = lysineCount + shiftNum;
                    expression = expression + "[Lysine Count] = " + K_to_Add;
                    if (shiftNum < (Convert.ToInt32(nUD_Missed_Ks.Value)))
                    {
                        expression = expression + " or ";
                    }
                    else
                    {
                        expression = expression + ") and ";
                    }
                }

                expression = expression + "[Light Retention Time] >= " + (retTime - Convert.ToDouble(nUD_RetTimeToleranace.Value)) +
                    " and [Light Retention Time] <= " + (retTime + Convert.ToDouble(nUD_RetTimeToleranace.Value));

                DataRow[] aggTheseRows = GlobalData.acceptableNeuCodeLightProteoforms.Select(expression);
                
                double aggIntSum = 0;

                foreach (DataRow row in aggTheseRows)
                {
                    //aggMassSum = aggMassSum + double.Parse(row["Light Mass Corrected"].ToString());
                    aggIntSum = aggIntSum + double.Parse(row["Light Intensity"].ToString());
                    //aggRTSum + double.Parse(row["Light Retention Time"].ToString());
                }

                double aggMass = 0;
                //double aggInt = 0;
                double aggRT = 0;

                foreach (DataRow row in aggTheseRows)
                {
                    double massShift = Math.Round((mass - double.Parse(row["Light Mass Corrected"].ToString())), 0) * 1.0015;
                    aggMass = aggMass + (double.Parse(row["Light Mass Corrected"].ToString())+massShift)*(double.Parse(row["Light Intensity"].ToString())/aggIntSum);
                    //aggInt = aggInt + double.Parse(row["Light Intensity"].ToString()) * (double.Parse(row["Light Intensity"].ToString()) / aggIntSum);
                    aggRT = aggRT + double.Parse(row["Light Retention Time"].ToString()) * (double.Parse(row["Light Intensity"].ToString()) / aggIntSum);
                }

                foreach (DataRow row in aggTheseRows)
                {
                    //MessageBox.Show("supposed to be updating table");
                    row["Aggregated Mass"] = aggMass;
                    row["Aggregated Intensity"] = aggIntSum;
                    row["Aggregated Retention Time"] = aggRT;
                }

                GlobalData.aggregatedProteoforms.Rows.Add(aggMass, aggIntSum, aggRT, lysineCount);
                //MessageBox.Show("agg proteoforms row added I think");
            }
            FillAggregatesTable();
        }

        

        private void dgv_AggregatedProteoforms_CellContentClick(object sender, DataGridViewCellEventArgs e)
        {
            double mass;

            if(e.RowIndex >= 0)
            {
                DataGridViewRow row = this.dgv_AggregatedProteoforms.Rows[e.RowIndex];
                mass = Convert.ToDouble(row.Cells["Aggregated Mass"].Value.ToString());

                DataTable dtClone = GlobalData.acceptableNeuCodeLightProteoforms.Clone();
                foreach (DataRow aRow in GlobalData.acceptableNeuCodeLightProteoforms.Select(
                    "[Aggregated Mass] > " + (mass-.001) + " and [Aggregated Mass] < " + (mass + .001)))
                {
                    dtClone.ImportRow(aRow);
                }

                //Round decimals before displaying
                string[] rt_column_names = new string[] { "Aggregated Retention Time", "Light Retention Time" };
                string[] intensity_column_names = new string[] { "Aggregated Intensity", "Light Intensity" };
                string[] mass_column_names = new string[] { "Aggregated Mass", "Light Mass Corrected", "Light Mass" };
                dataTableHandler.DisplayWithRoundedDoubles(dgv_AcceptNeuCdLtProteoforms, dtClone,
                    rt_column_names, intensity_column_names, new string[] { }, mass_column_names);
            }
        }

        private void nUP_mass_tolerance_ValueChanged(object sender, EventArgs e)
        {
            if (GlobalData.acceptableNeuCodeLightProteoforms.Rows.Count > 0)
            {
                ZeroAggregateMasses();
                AggregateNeuCodeLightProteoforms();
            }
        }

        private void nUD_RetTimeToleranace_ValueChanged(object sender, EventArgs e)
        {
            if (GlobalData.acceptableNeuCodeLightProteoforms.Rows.Count > 0)
            {
                ZeroAggregateMasses();
                AggregateNeuCodeLightProteoforms();
            }
        }

        private void nUD_Missed_Monos_ValueChanged(object sender, EventArgs e)
        {
            if (GlobalData.acceptableNeuCodeLightProteoforms.Rows.Count > 0)
            {
                ZeroAggregateMasses();
                AggregateNeuCodeLightProteoforms();
            }
        }

        private void nUD_Missed_Ks_ValueChanged(object sender, EventArgs e)
        {
            if (GlobalData.acceptableNeuCodeLightProteoforms.Rows.Count > 0)
            {
                ZeroAggregateMasses();
                AggregateNeuCodeLightProteoforms();
            }
        }
    }
}
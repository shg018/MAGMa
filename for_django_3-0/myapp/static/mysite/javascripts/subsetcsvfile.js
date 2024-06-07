var sliderFC = document.getElementById("fcRange");
// var outputFC = document.getElementById("fcValue");
// outputFC.innerHTML = sliderFC.value; // Display the default FC slider value

// // Update the FC slider value (each time you drag the slider handle)
// sliderFC.oninput = function () {
//     outputFC.innerHTML = this.value;
// }

var sliderPval = document.getElementById("pvalRange");
// var outputPval = document.getElementById("pvalValue");
// outputPval.innerHTML = sliderPval.value; // Display the default FC slider value

// // Update the Pvalue slider value (each time you drag the slider handle)
// sliderPval.oninput = function () {
//     outputPval.innerHTML = this.value;
// }

var sliderPSM = document.getElementById("psmRange");
// var outputPSM = document.getElementById("psmValue");
// outputPSM.innerHTML = sliderPSM.value; // Display the default FC slider value

// // Update the PSM slider value (each time you drag the slider handle)
// sliderPSM.oninput = function () {
//     outputPSM.innerHTML = this.value;
// }




document.addEventListener("DOMContentLoaded", function () {
    // const slider1 = document.getElementById("slider1");
    // const slider2 = document.getElementById("slider2");

    // const downloadButton = document.getElementById("download-subset-file-button");
    var filename_volcano_plot = document.getElementById('volcanoplot-file').innerHTML;
    var filename_for_save = filename_volcano_plot.replace("Results/", "");
    filename_for_save = filename_for_save.replace(".csv", "");
    var tmpstr = filename_for_save.split("/");
    filename_for_save = tmpstr[3];
    // alert(filename_for_save);
    // var filename_for_save = keeptmpstr;
    filename_volcano_plot = filename_volcano_plot.replace("Results/", "static/");

    let fcval = +Math.log2(sliderFC.value);
    let pvalcutoff = sliderPval.value;
    let psmcutoff = sliderPSM.value;
    alert(fcval);


    const downloadButton = document.getElementById("download-subset-file-button");

    function updateSliderValues() {
        fcval = +Math.log2(sliderFC.value);
        pvalcutoff = sliderPval.value;
        psmcutoff = sliderPSM.value;
    }

    // var sliderFC = document.getElementById("fcRange");
    // var sliderPval = document.getElementById("pvalRange");
    // var sliderPSM = document.getElementById("psmRange");

    // var filename_volcano_plot = document.getElementById('volcanoplot-file').innerHTML;
    // var filename_for_save = filename_volcano_plot.replace("Results/", "");
    // filename_for_save = filename_for_save.replace(".csv", "");
    // var tmpstr = filename_for_save.split("/");
    // filename_for_save = tmpstr[3];
    // alert(filename_for_save);
    // // var filename_for_save = keeptmpstr;
    // filename_volcano_plot = filename_volcano_plot.replace("Results/", "static/");

    // var fcval = +sliderFCval;
    // alert(fcval);
    // var pvalcutoff = +sliderPval;
    // var psmcutoff = +sliderPSM;

    // Function to generate CSV content based on slider values
    function generateCSVContent(data) {
        // Perform filtering based on slider values
        const filteredData = data.filter(function (d) {
            return d.log2FC >= fcval && d.adjpval <= +pvalcutoff && d.numpsms >= +psmcutoff;
        });

        // Convert filtered data back to CSV format
        const csvContent = d3.csvFormat(filteredData);
        return csvContent;
    }
    

    // Function to read CSV file and subset data
    function readAndSubsetCSV() {
        // fcval,pvalcutoff,psmcutoff
        d3.csv(filename_volcano_plot, function (data) {
            

            data.forEach(function (d) {
                d.pval = +d.pval;
                d.adjpval = +d.adjpval;
                d.neglogpval = -1.0 * Math.log10(d.pval);
                d.log2FC = +d.log2FC;
                d.genesymbol = d['Gene Symbol'];
                d.Protein = d['Protein'];
                d.numpsms = +d['# PSMs'];
            });

            // // Perform filtering based on slider values
            // const filteredData = data.filter(function (d) {
            //     return d.log2FC >= fcval && d.adjpval <= pvalcutoff && d.numpsms >= psmcutoff;
            // });
            // Convert filtered data back to CSV format
            // const csvContent = d3.csvFormat(filteredData);
            const csvContent = generateCSVContent(data);

            // Download the subset CSV file
            downloadCSV(filename_for_save, csvContent);
        });
        // .catch(function (error) {
        //     console.error("Error reading CSV:", error);
        // });
    }

    // Function to download CSV file
    function downloadCSV(filename_for_save, csvContent) {
        const blob = new Blob([csvContent], { type: "text/csv" });
        const url = window.URL.createObjectURL(blob);
        const a = document.createElement("a");
        a.href = url;
        a.download = filename_for_save + "_subset_enriched.csv";
        document.body.appendChild(a);
        a.click();
        window.URL.revokeObjectURL(url);
        document.body.removeChild(a);
    }

    // Event listener for download button click
    downloadButton.addEventListener("click", function () {
        // sliderFC.addEventListener("change", function () {
        //     fcval = this.value;
        //     alert(fcval);
        // });
    
        // sliderPval.addEventListener("change", function () {
        //     pvalcutoff = this.value;
        // });
    
        // sliderPSM.addEventListener("change", function () {
        //     psmcutoff = this.value;
        // });
        // alert("Download button");
        // fcval,pvalcutoff,psmcutoff
        readAndSubsetCSV();
    });

    // Listen to the FC slider?
    sliderFC.addEventListener("input", function () {
        updateSliderValues();
    });

    sliderPval.addEventListener("input", function () {
        updateSliderValues();
    });

    sliderPSM.addEventListener("input", function () {
        updateSliderValues();
    });
});
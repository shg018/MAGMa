var margin = { top: 50, right: 30, bottom: 30, left: 60 },
    widthvolc = 860 - margin.left - margin.right,
    heightvolc = 600 - margin.top - margin.bottom;



// append the svg object to the body of the page
var svg = d3.select("#my_dataviz")
    .append("svg")
    .attr("width", widthvolc + margin.left + margin.right)
    .attr("height", heightvolc + margin.top + margin.bottom)
    .append("g")
    .attr("transform",
        "translate(" + margin.left + "," + margin.top + ")");


var sliderFC = document.getElementById("fcRange");
var outputFC = document.getElementById("fcValue");
outputFC.innerHTML = sliderFC.value; // Display the default FC slider value

// Update the FC slider value (each time you drag the slider handle)
sliderFC.oninput = function () {
    outputFC.innerHTML = this.value;
}

var sliderPval = document.getElementById("pvalRange");
var outputPval = document.getElementById("pvalValue");
outputPval.innerHTML = sliderPval.value; // Display the default FC slider value

// Update the Pvalue slider value (each time you drag the slider handle)
sliderPval.oninput = function () {
    outputPval.innerHTML = this.value;
}

var sliderPSM = document.getElementById("psmRange");
var outputPSM = document.getElementById("psmValue");
outputPSM.innerHTML = sliderPSM.value; // Display the default FC slider value

// Update the PSM slider value (each time you drag the slider handle)
sliderPSM.oninput = function () {
    outputPSM.innerHTML = this.value;
}

var sliderDotsize = document.getElementById("dotsizeRange");
var outputDotsize = document.getElementById("dotsizeValue");
outputDotsize.innerHTML = sliderDotsize.value; // Display the default FC slider value

// Update the dot size slider value (each time you drag the slider handle)
sliderDotsize.oninput = function () {
    outputDotsize.innerHTML = this.value;
}

//Read the data
var filename_volcano_plot = document.getElementById('volcanoplot-file').innerHTML;
var filename_for_save = filename_volcano_plot.replace("Results/", "");
filename_for_save = filename_for_save.replace(".csv", "");
var tmpstr = filename_for_save.split("/");
var tmpstr2 = tmpstr[3].split("_FC_");
filename_for_save = tmpstr2[0];
filename_for_save = filename_for_save.replace("-vs-", " / ");
//TODO - static file not found and app dies as a result
filename_volcano_plot = filename_volcano_plot.replace("Results/", "static/");

d3.csv(filename_volcano_plot, function (data) {

    data.forEach(function (d) {
        d.pval = +d.pval;
        d.adjpval = +d.adjpval;
        d.neglogpval = -1.0 * Math.log10(d.pval);
        d.log2FC = +d.log2FC;
        d.genesymbol = d['Gene Symbol'];
        d.Protein = d['Protein'];
        d.numpsms = +d['# PSMs'];
        // d.genesymbol = parseDate(d['Gene Symbol']);

    });
    // alert(data);

    // #TODO filter data by # PSMs
    // data.filter(function(d){ return  (d.name == "toto" || d.name == "tutu") })

    var height_yaxis = d3.max(data, function (d) { return d.neglogpval; });
    var pos_width_xaxis = d3.max(data, function (d) { return d.log2FC; });
    var neg_width_xaxis = d3.min(data, function (d) { return d.log2FC; });
    var x_range = d3.max([pos_width_xaxis, -1.0 * neg_width_xaxis]);
    // alert(height_yaxis);
    // alert(pos_width_xaxis);
    // alert(x_range);

    // Add X axis
    var x = d3.scaleLinear()
        .domain([-1.0 * x_range, x_range])
        // .domain([0, 1000])
        .range([0, widthvolc]);
    svg.append("g")
        .attr("transform", "translate(0," + heightvolc + ")")
        .call(d3.axisBottom(x));


    // Add Y axis - shift to left
    var y = d3.scaleLinear()
        .domain([0, height_yaxis])
        // .domain([0,800])
        .range([heightvolc, x_range]);
    svg.append("g")
        .attr("transform", "translate(" + 1 + ",0)")
        .call(d3.axisLeft(y));

    // Append axes to the chart
    svg.append("g")
        .attr("transform", `translate(0, ${heightvolc})`)
        .call(x);

    svg.append("g")
        .call(y);

    // Create title on the volcano plot
    svg.append("text")
        .attr("transform", `translate(${widthvolc / 2}, ${5})`)
        .attr("class", "titleText")
        .style("text-anchor", "middle")
        .text(filename_for_save)
        .style("font-size","20px");

    // Create axes labels
    svg.append("text")
        .attr("transform", "rotate(-90)")
        .attr("y", 0 - margin.left + 20)
        .attr("x", 0 - (heightvolc / 2))
        .attr("dy", "1em")
        .attr("class", "axisText")
        .style("text-anchor", "middle")
        .text("-log10(P-value)");

    svg.append("text")
        .attr("transform", `translate(${widthvolc / 2}, ${heightvolc + 25})`)
        .attr("class", "axisText")
        .style("text-anchor", "middle")
        .text("log2FC");


    svg.append("ygridline")
        .attr("x1", 10)  //<<== change your code here
        .attr("y1", 0)
        .attr("x2", 50)  //<<== and here
        .attr("y2", 50)
        // .attr("x1", x(widthvolc/2))  //<<== change your code here
        // .attr("y1", 0)
        // .attr("x2", x(widthvolc/2))  //<<== and here
        // .attr("y2", heightvolc)
        .style("stroke-width", 10)
        .style("stroke", "red")
    //.style("fill", "none");

    // Add grid lines AND FDR thresholds (10%, 5%, 1%)


    // alert(height_yaxis);
    // alert(outputFC.innerHTML);

    var FCval = Math.log2(outputFC.innerHTML);
    var pvalcutoff = outputPval.innerHTML;//-1.0*Math.log10(outputPval.innerHTML);
    var psmcutoff = outputPSM.innerHTML;
    var dotsize = outputDotsize.innerHTML;
    // alert(dotsize);

    // alert(FCval);

    var tooltip = d3.select("#my_dataviz")
        .append("g")
        // .attr("class","custom-tooltip-tooltiptext")
        .style("position", "absolute")
        .style("z-index", "500")
        .style("visibility", "hidden")
        .style("background", "#969696")
        .style("color", "#fff")
        .style("border", "3px solid #F15A24")
        .style("padding", "5px 5px");
    // .text("a simple tooltip");


    function hexcodeScatter(d, FCval, pvalcutoff) {
        if ((d.log2FC >= FCval) && (d.adjpval <= pvalcutoff)) { return "#31a354"; }
        else if ((d.log2FC >= FCval) && (d.adjpval > pvalcutoff)) { return "#a1d99b"; }
        else if (((d.log2FC < FCval) && (d.log2FC > -1.0 * FCval)) && (d.adjpval <= pvalcutoff)) { return "#636363"; }
        else if (((d.log2FC < FCval) && (d.log2FC > -1.0 * FCval)) && (d.adjpval > pvalcutoff)) { return "#bdbdbd"; }
        else if ((d.log2FC <= -1.0 * FCval) && (d.adjpval <= pvalcutoff)) { return "#de2d26"; }
        else { return "#fc9272"; }

    }

    var startData = data.filter(function (d) { return d.numpsms >= psmcutoff })
    var startDotsize = dotsize
    var uniprotInput = document.getElementById("protein-list-input").value.split("@");
    // Add dots for scatter plot
    var scatter = svg.selectAll("circle")
        .data(startData)
        //.filter(function(d) {return (d.numpsms > psmcutoff) })
        .enter()
        .append("circle")
        .attr("cx", function (d) { return x(d.log2FC); })
        .attr("cy", function (d) { return y(d.neglogpval); })
        .attr("r", startDotsize)
        .style("fill", function (d) { return hexcodeScatter(d, FCval, pvalcutoff); })
        //.text(function(d) { return d.genesymbol; })
        .on("mouseover", function (d) { tooltip.style("border", "3px solid " + hexcodeScatter(d, FCval, pvalcutoff)); tooltip.html((d.genesymbol + "<br>" + "Protein: " + d.Protein + "<br>" + "# PSMs: " + d.numpsms + "<br>log2FC: " + (d.log2FC).toFixed(2) + "<br>P-value: " + Number.parseFloat(d.pval).toExponential(2) + "<br>Adj.P-value: " + Number.parseFloat(d.adjpval).toExponential(2))); tooltip.style("top", (event.pageY - d.neglogpval + 3) + "px"); tooltip.style("left", (event.pageX - d.log2FC + 3) + "px"); return tooltip.style("visibility", "visible"); })
        .on("mouseout", function () { return tooltip.style("visibility", "hidden"); });

    // svg.append("circletext")

    svg.selectAll("g")
        .data(startData)
        .enter()
        .append("text")
        .text(function (d) { return ""; })
        .attr("x", function (d) { return x(d.log2FC) - 2; }) // Adjust the positioning as needed
        .attr("y", function (d) { return y(d.neglogpval) + 2; }) // Adjust the positioning as needed
        .attr("class", "scatterlabel");  // Add a class for styling

    function updatePointSize(newSize) {
        // Update the point size
        svg.selectAll("circle").attr("r", newSize);
    }

    // A function that update the chart when slider is moved?
    function updateChart(FCval, pvalcutoff, psmcutoff, dotsize, uniprotInput) {
        // // update the chart
        var newData = data.filter(function (d) { return d.numpsms >= psmcutoff })
        update = svg.selectAll("circle")
            .data(newData)
        // Remove some of the points filtered as a result
        update.exit().remove();
        // Add the new points as a result of the filter
        update.enter()
            .append("circle")
            .attr("cx", function (d) { return x(d.log2FC); })
            .attr("cy", function (d) { return y(d.neglogpval); })
            .attr("r", dotsize)
            .style("fill", function (d) { return hexcodeScatter(d, FCval, pvalcutoff); })
            .on("mouseover", function (d) { tooltip.style("border", "3px solid " + hexcodeScatter(d, FCval, pvalcutoff)); tooltip.html((d.genesymbol + "<br>" + "Protein: " + d.Protein + "<br>" + "# PSMs: " + d.numpsms + "<br>log2FC: " + (d.log2FC).toFixed(2) + "<br>P-value: " + Number.parseFloat(d.pval).toExponential(2) + "<br>Adj.P-value: " + Number.parseFloat(d.adjpval).toExponential(2))); tooltip.style("top", (event.pageY - d.neglogpval + 3) + "px"); tooltip.style("left", (event.pageX - d.log2FC + 3) + "px"); return tooltip.style("visibility", "visible"); })
            .on("mouseout", function () { return tooltip.style("visibility", "hidden"); });

        // Keep the dots that should be retained even with the filter and update their values
        update.attr("cx", function (d) { return x(d.log2FC); })
            .attr("cy", function (d) { return y(d.neglogpval); })
            .attr("r", dotsize)
            .style("fill", function (d) { return hexcodeScatter(d, FCval, pvalcutoff); })
            .on("mouseover", function (d) { tooltip.style("border", "3px solid " + hexcodeScatter(d, FCval, pvalcutoff)); tooltip.html((d.genesymbol + "<br>" + "Protein: " + d.Protein + "<br>" + "# PSMs: " + d.numpsms + "<br>log2FC: " + (d.log2FC).toFixed(2) + "<br>P-value: " + Number.parseFloat(d.pval).toExponential(2) + "<br>Adj.P-value: " + Number.parseFloat(d.adjpval).toExponential(2))); tooltip.style("top", (event.pageY - d.neglogpval + 3) + "px"); tooltip.style("left", (event.pageX - d.log2FC + 3) + "px"); return tooltip.style("visibility", "visible"); })
            .on("mouseout", function () { return tooltip.style("visibility", "hidden"); });

        updatePlotLabels(newData, uniprotInput);


    }

    function updatePlotLabels(newData, labels) {
        // var circleLabels = svg.selectAll("g").data(newData);
        // circleLabels.select("scatterlabel").remove();

        for (let i = 0; i < labels.length; i++) {
            // P29590;Q9BQE3
            // var filteredData = newData.filter(function(d){ return d.Protein = labels[i]})
            // svg.selectAll("g").data(newData).enter().text(function (d) { if(d.Protein==labels[i]){return d.genesymbol;}else{return "";} })

            svg.selectAll("g")
                .data(newData)//.selectAll("text").text(function (d) { if(d.Protein==labels[i]){return d.genesymbol;}else{return "";} })
                .enter()
                .append("text")
                .text(function (d) { if (d.Protein == labels[i]) { return d.genesymbol; } else { return ""; } })
                .attr("x", function (d) { return x(d.log2FC) - 2; }) // Adjust the positioning as needed
                .attr("y", function (d) { return y(d.neglogpval) + 2; }) // Adjust the positioning as needed
                .attr("class", "label");  // Add a class for styling


            // circleLabels.enter()
            // circleLabels.enter().each(function (d) {
            //   // Remove existing label if any
            //   // d3.select(this).select("text").remove();
            //   // Add labels only to circles that have labels provided
            //   if (d.Protein==labels[i]) {
            //       d3.select(this)
            //       .append("text")
            //       .attr("x", function (d) { return x(d.log2FC)-2; })
            //       .attr("y", function (d) { return y(d.neglogpval)+2; }) // Adjust the positioning as needed
            //       .text(d.genesymbol)
            //       .attr("class", "scatterlabel");
            //   }
            // });
        }
    }

    // Listen to the FC slider?
    d3.select("#fcRange").on("change", function (d) {
        selectedValue = this.value
        //console.log(selectedValue)
        FCval = Math.log2(selectedValue)
        updateChart(FCval, pvalcutoff, psmcutoff, dotsize, uniprotInput)
    });

    // Listen to the Pvalue slider?
    d3.select("#pvalRange").on("change", function (d) {
        selectedValue = this.value
        //console.log(selectedValue)
        pvalcutoff = selectedValue//-1.0*Math.log10(selectedValue)
        updateChart(FCval, pvalcutoff, psmcutoff, dotsize, uniprotInput)
    });

    // Listen to the PSM slider?
    d3.select("#psmRange").on("change", function (d) {
        selectedValue = this.value
        //console.log(selectedValue)
        psmcutoff = selectedValue
        updateChart(FCval, pvalcutoff, psmcutoff, dotsize, uniprotInput)
    });

    // Listen to the Dot size slider?
    d3.select("#dotsizeRange").on("change", function (d) {
        selectedValue = this.value
        dotsize = selectedValue
        updatePointSize(dotsize)
    });

    // Listen to the Update button to assign labels to scatter plot
    d3.select("#update-button").on("click", function (d) {
        newData = data.filter(function (d) { return d.numpsms >= psmcutoff })
        uniprotInput = document.getElementById("protein-list-input").value.split("@");
        // alert(uniprotInput)
        updatePlotLabels(newData, uniprotInput)
    });


    // element = d3.select(".chart").style("font", "10px sans-serif");
});

// d3.select(".chart").style("font", "10px sans-serif");
// document.getElementById("download-button").disabled = false;

// function hyperbolicRightCurve(X,c,x0plot){
//     return X.map(function(x) {
//         return [x, c/(i-2*x0plot)];
//     });
// }

// console.log(element.outerHTML);
function SVG2PNG(objsvg, callback) {
    var canvas = document.createElement('canvas'); // Create a Canvas element.
    var ctx = canvas.getContext('2d'); // For Canvas returns 2D graphic.
    var data = objsvg.outerHTML;//new XMLSerializer().serializeToString(svg.node());//svg.outerHTML; // Get SVG element as HTML code.
    canvg(canvas, data); // Render SVG on Canvas.
    callback(canvas); // Execute callback function.
}
function generateLink(fileName, data) {
    var link = document.createElement('a');
    link.download = fileName;
    link.href = data;
    return link;
}
document.getElementById('download-button').onclick = function () { // Bind click event on download button.
    var element = document.getElementById('my_dataviz').getElementsByTagName('svg')[0]; // Get SVG element.

    SVG2PNG(element, function (canvas) { // Arguments: SVG element, callback function.
        var base64 = canvas.toDataURL("image/png"); // toDataURL return DataURI as Base64 format.
        generateLink('volcano_plot.png', base64).click(); // Trigger the Link is made by Link Generator and download.
    });
}
// document.getElementById('download-button').addEventListener('click', function () {
//   // Get the SVG element
//   var svgElement = document.getElementById('my_dataviz');//.getElementsByTagName('svg')[0];
//   var canvas = document.createElement('canvas'); 
//   // Create a canvas from the SVG
//   html2canvas(svgElement).then(function (canvas) {
//     // Convert the canvas to a data URL
//     var dataURL = canvas.toDataURL('image/jpeg');

//     // Create a PDF document
//     var pdf = new jsPDF();

//     // window.jsPDF = window.jspdf.jsPDF;

//     // Add an image to the PDF
//     pdf.addImage(dataURL, 'JPEG');//, 0, 0, svgElement.width.baseVal.value, svgElement.height.baseVal.value);

//     // Save the PDF
//     pdf.save('volcano_plot.pdf');
//   });
// })

function showTooltipRight(evt, text) {
    let tooltip = document.getElementById("tooltip");
    tooltip.innerHTML = text;
    tooltip.style.display = "block";
    tooltip.style.left = evt.pageX - 150 + 'px';
    tooltip.style.top = evt.pageY + 10 + 'px';
}

function showTooltipLeft(evt, text) {
    let tooltip = document.getElementById("tooltip");
    tooltip.innerHTML = text;
    tooltip.style.display = "block";
    tooltip.style.left = evt.pageX + 150 + 'px';
    tooltip.style.top = evt.pageY + 10 + 'px';
}

function hideTooltip() {
    var tooltip = document.getElementById("tooltip");
    tooltip.style.display = "none";
}
// }



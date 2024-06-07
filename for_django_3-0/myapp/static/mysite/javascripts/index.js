// Set smaller font size if viewport is small (e.g. iPhone 5)
if (window.innerWidth < 375) {
    $(".footer-text").css("font-size", "8pt")
}

// Arc data
let data = [
    {name: "nsp1", value: 1},
    {name: "nsp2", value: 1},
    {name: "nsp3", value: 1},
    {name: "nsp4", value: 1},
    {name: "nsp5", value: 1},
    {name: "nsp6", value: 1},
    {name: "nsp7", value: 1},
    {name: "nsp8", value: 1},
    {name: "nsp9", value: 1},
    {name: "nsp10", value: 1},
    {name: "nsp11", value: 1},
    {name: "nsp12", value: 1},
    {name: "nsp13", value: 1},
    {name: "nsp14", value: 1},
    {name: "nsp15", value: 1},
    {name: "nsp16", value: 1},
    {name: "Spike", value: 1},
    {name: "orf3a", value: 1},
    {name: "orf3b", value: 1},
    {name: "E", value: 1},
    {name: "M", value: 1},
    {name: "orf6", value: 1},
    {name: "orf7a", value: 1},
    {name: "orf7b", value: 1},
    {name: "orf8", value: 1},
    {name: "N", value: 1},
    {name: "orf9b", value: 1},
    {name: "orf9c", value: 1},
    {name: "orf10", value: 1},
    {name: "", value: 1}, // Blank
];
for (let i = 0; i < data.length; i++) {
    if (data[i].name !== "") {
        let o = new Option(data[i].name.toUpperCase(), "COVID19" + data[i].name);
        $(o).html(data[i].name.toUpperCase());
        $("#prot1").append(o);
    }
}

// Global variable
let currChosen = undefined;

// Populate the svg container.
function drawShapes(id) {
    // Drawing donut chart.
    let color = d3.scaleOrdinal()
        .domain(data.map(d => d.name))
        .range(d3.quantize(t => d3.interpolateGnBu(t * 0.6 + 0.2), data.length - 1));
    let svgWidth = $("#" + id).width();
    let svgHeight = Math.min($("#" + id).width(), $("main").innerHeight());
    let svgContainer = d3.select("#" + id).append("svg").attr('width', svgWidth).attr('height', svgHeight);
    let pie = d3.pie().padAngle(0.005).sort(null).value(d => d.value);
    let arcs = pie(data);
    let arc = d3.arc()
        .innerRadius(svgHeight / 2 * 0.6)
        .outerRadius(svgHeight / 2 * 0.8)
        .startAngle(function(d) {return d.startAngle + 6 * Math.PI / 180;})
        .endAngle(function(d) {return d.endAngle + 6 * Math.PI / 180;});

    // Adding arcs.
    svgContainer.selectAll("path")
        .data(arcs.filter(function (d) {
            if (d.data.name === "") return false;
            return true;
        }))
        .join("path")
          .attr("fill", function(d) {return color(d.data.name);})
          .attr("id", d => `${d.data.name}-arc`)
          .attr("class", function(d) {if (d.data.name === "") return ""; else return "prot-arc";})
          .attr("d", arc)
          .attr("transform", "translate(" + (svgWidth / 2).toString() + ","  + (svgHeight / 2).toString() + ")")
          .style("opacity", 1)
        .append("title")
          .text(d => `${d.data.name}`);

    let circle = svgContainer.append("circle")
        .attr("cx", svgWidth / 2)
        .attr("cy", svgHeight / 2)
        .attr("r", svgHeight / 2 * 0.45)
        .attr("class", "svg-circle")
        .style("fill", "#eaf5f5");

    // Adding text.
    svgContainer.append("g")
            .attr("font-size", "0.8rem")
            .attr("text-anchor", "middle")
            .attr("transform", "translate(" + (svgWidth / 2).toString() + "," + (svgHeight / 2).toString() + ")")
        .selectAll("text")
        .data(arcs)
        .join("text")
          .attr("transform", d => `translate(${arc.centroid(d)})`)
          .call(text => text.append("tspan")
              .attr("pointer-events", "none")
              .attr("y", "0")
              .text(d => d.data.name.toUpperCase()));

    // Adding polylines.
    let borderArc = d3.arc()
        .innerRadius(svgHeight / 2 * 0.85)
        .outerRadius(svgHeight / 2 * 0.85)
        .startAngle(function(d) {return d.startAngle + 6 * Math.PI / 180;})
        .endAngle(function(d) {return d.endAngle + 6 * Math.PI / 180;});

    let outerArc = d3.arc()
        .innerRadius(svgHeight / 2 * 0.95)
        .outerRadius(svgHeight / 2 * 0.95)
        .startAngle(function(d) {return d.startAngle + 6 * Math.PI / 180;})
        .endAngle(function(d) {return d.endAngle + 6 * Math.PI / 180;});

    svgContainer.selectAll('allPolylines')
        .data(arcs.filter(function (d) {
            if (d.data.name === "") return false;
            return noStructure.indexOf(d.data.name) === -1;
        })).enter()
        .append('polyline')
        .attr("class", function(d) {if (d.data.name === "") return ""; else return "prot-polyline";})
        .attr("stroke", "#cccccc")
        .style("fill", "none")
        .attr("stroke-width", 2)
        .attr("transform", "translate(" + (svgWidth / 2).toString() + ","  + (svgHeight / 2).toString() + ")")
        .attr("id", d => `${d.data.name}-polyline`)
        .attr("hidden", true)
        .attr('points', function(d) {
            let posA = borderArc.centroid(d); // line insertion in the slice
            let posB = outerArc.centroid(d); // line break: we use the other arc generator that has been built only for that
            let posC = outerArc.centroid(d); // Label position = almost the same as posB
            let midangle = d.startAngle + (d.endAngle - d.startAngle) / 2; // we need the angle to see if the X position will be at the extreme right or extreme left
            posC[0] = 1.1 * svgHeight / 2 * 0.95 * (midangle < Math.PI ? 1 : -1); // multiply by 1 or -1 to put it on the right or on the left
            return [posA, posB, posC];
        });

    // Adding images.
    svgContainer.selectAll("image")
        .data(arcs.filter(function (d) {
            if (d.data.name === "") return false;
            return noStructure.indexOf(d.data.name) === -1;
        })).enter()
        .append("image")
        .attr("xlink:href", function (d) {return staticPrefix + "protein_images/" + d.data.name.toUpperCase() + "_Sphere_Transparent.png";})
        .attr("height", 200).attr("width", 200)
        .attr("id", d => `${d.data.name}-image`)
        .attr("hidden", true)
        .attr('transform', function(d) {
            let pos = outerArc.centroid(d);
            let midangle = d.startAngle + (d.endAngle - d.startAngle) / 2;
            if (midangle < Math.PI) pos[0] = 1.1 * svgHeight / 2 * 0.95 + (svgWidth / 2);
            else pos[0] = -1.1 * svgHeight / 2 * 0.95 + (svgWidth / 2) - 200;
            pos[1] = pos[1] + (svgHeight / 2) - 100;
            if (pos[1] <= 0) pos[1] = pos[1] + 60;
            else if (pos[1] + 200 >= svgHeight) pos[1] = pos[1] - 60;
            return 'translate(' + pos + ')';
        });

    // Handling mouse events.
    svgContainer.selectAll("path")
        .on("mouseover", function(d, i, j) {
            d3.select(this)
                .transition()
                .attr('d', d3.arc()
                    .innerRadius(svgHeight / 2 * 0.6)
                    .outerRadius(svgHeight / 2 * 0.85)
                    .startAngle(function(d) {return d.startAngle + 6 * Math.PI / 180;})
                    .endAngle(function(d) {return d.endAngle + 6 * Math.PI / 180;})
                )
                .style("fill", function(d) {
                    if (d.data.name === "") return "#f5f5f5";
                    else return d3.color(d3.select(this).attr("fill")).darker(1)
                });
            $(`#${this.id.split('-')[0]}-polyline`).removeAttr("hidden");
            $(`#${this.id.split('-')[0]}-image`).removeAttr("hidden");
        })
        .on("mouseout", function(d, i, j) {
            if (! $(`#${this.id.split('-')[0]}-arc`).hasClass("arc-clicked")) {
                d3.select(this)
                    .transition()
                    .attr('d', d3.arc()
                        .innerRadius(svgHeight / 2 * 0.6)
                        .outerRadius(svgHeight / 2 * 0.8)
                        .startAngle(function(d) {return d.startAngle + 6 * Math.PI / 180;})
                        .endAngle(function(d) {return d.endAngle + 6 * Math.PI / 180;})
                    )
                    .style("fill", function (d) {
                        if (d.data.name === "") return "#f5f5f5";
                        else return color(d.data.name);
                    });
                $(`#${this.id.split('-')[0]}-polyline`).attr("hidden", true);
                $(`#${this.id.split('-')[0]}-image`).attr("hidden", true);
            }
        })
        .on("click", function(d, i, j){
            if (this.id.split('-')[0] === currChosen) return;
            if (currChosen !== undefined) {
                let oldArcDiv = $(`#${currChosen}-arc`);
                d3.select(oldArcDiv[0])
                    .transition()
                    .attr('d', d3.arc()
                        .innerRadius(svgHeight / 2 * 0.6)
                        .outerRadius(svgHeight / 2 * 0.8)
                        .startAngle(function(d) {return d.startAngle + 6 * Math.PI / 180;})
                        .endAngle(function(d) {return d.endAngle + 6 * Math.PI / 180;})
                    )
                    .style("fill", function (d) {
                        if (d.data.name === "") return "#f5f5f5";
                        else return color(d.data.name);
                    });
                $(`#${currChosen}-polyline`).attr("hidden", true);
                $(`#${currChosen}-image`).attr("hidden", true);
                oldArcDiv.removeClass("arc-clicked");
            }
            let arcDiv = $(`#${this.id.split('-')[0]}-arc`);
            if (arcDiv.hasClass("arc-clicked")) {
                arcDiv.removeClass("arc-clicked");
                currChosen = undefined;
                $("#prot1").val("--");
                populateInteractor();
            } else {
                arcDiv.addClass("arc-clicked");
                currChosen = this.id.split('-')[0];
                $("#prot1").val("COVID19" + currChosen);
                populateInteractor();
            }
        });
}

// Dropdown menu.
function populateInteractor() {
    $("#prot2").empty();
    let selectedViral = $("#prot1").val();
    let disabledAttr = $("#prot2").attr("disabled");
    if (selectedViral === "--" || interactionByViralProt[selectedViral.slice(7)] === undefined) {
        if (typeof disabledAttr === typeof undefined || disabledAttr === false) $("#prot2").attr("disabled", true);
        return;
    }
    for (let i = 0; i < interactionByViralProt[selectedViral.slice(7)].length; i++) {
        if (interactionByViralProt[selectedViral.slice(7)][i] !== "") {
            let o = new Option(interactionByViralProt[selectedViral.slice(7)][i], interactionByViralProt[selectedViral.slice(7)][i]);
            $(o).html(humanUniprotToGene[interactionByViralProt[selectedViral.slice(7)][i]]);
            $("#prot2").append(o);
        }
    }
    if (typeof disabledAttr !== typeof undefined && disabledAttr !== false) $("#prot2").removeAttr("disabled");
}

// Draw shapes and set form dimensions.
drawShapes("arc-svg");
let innerDiameter = Math.min($("#arc-svg").width(), $("main").innerHeight()) * 0.45;
$("#form-div")
    .css("height", innerDiameter / Math.sqrt(2))
    .css("width", innerDiameter / Math.sqrt(2))
    .css("margin-left", -innerDiameter / (Math.sqrt(2) * 2))
    .css("margin-top", -innerDiameter / (Math.sqrt(2) * 2));

// Handles window resize.
$(window).on('resize', function() {
    d3.select("#arc-svg > svg").remove();
    drawShapes("arc-svg");
    innerDiameter = Math.min($("#arc-svg").width(), $("main").innerHeight()) * 0.45;
    $("#form-div")
        .css("height", innerDiameter / Math.sqrt(2))
        .css("width", innerDiameter / Math.sqrt(2))
        .css("margin-left", -innerDiameter / (Math.sqrt(2) * 2))
        .css("margin-top", -innerDiameter / (Math.sqrt(2) * 2));
});

// Particle JS background.
particlesJS("particles-js", {
    "particles": {
        "number":{
            "value":80,
            "density":{"enable":true,"value_area":800},
        },
        "color":{"value":"#ff3b00"},
        "shape":{
            "type":"circle",
            "stroke":{"width":0,"color":"#000000"},
            "polygon":{"nb_sides":5},
            "image":{"src":"img/github.svg","width":100,"height":100},
        },
        "opacity":{
            "value":0.5,
            "random":false,
            "anim":{
                "enable":false,
                "speed":1,
                "opacity_min":0.1,
                "sync":false},
        },
        "size":{
            "value":3,
            "random":true,
            "anim":{
                "enable":false,
                "speed":40,
                "size_min":0.1,
                "sync":false
            },
        },
        "line_linked":{
            "enable":true,
            "distance":150,
            "color":"#e32828",
            "opacity":0.4,
            "width":1
        },
        "move":{
            "enable":true,
            "speed":2,
            "direction":"none",
            "random":false,
            "straight":false,
            "out_mode":"out",
            "bounce":false,
            "attract":{"enable":false,"rotateX":600,"rotateY":1200},
        },
    },
    "interactivity":{
        "detect_on":"canvas",
        "events":{
            "onhover":{"enable":false,"mode":"repulse"},
            "onclick":{"enable":true,"mode":"push"},
            "resize":true
        },
        "modes":{
            "grab":{"distance":400,"line_linked":{"opacity":1}},
            "bubble":{"distance":400,"size":40,"duration":2,"opacity":8,"speed":3},
            "repulse":{"distance":200,"duration":0.4},"push":{"particles_nb":4},
            "remove":{"particles_nb":2},
        },
    },
    "retina_detect":true});
update = function() {
    requestAnimationFrame(update);
};
requestAnimationFrame(update);
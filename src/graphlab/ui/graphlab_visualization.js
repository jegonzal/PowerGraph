google.load("jquery", "1.4.2");
google.load("jqueryui", "1.7.2");
google.load("visualization", "1", 
            {"packages":["corechart", "table", "gauge"]});


var domain_str = "http://localhost:8090/"
var update_interval = 2000;


// Start the rendering of the UI
google.setOnLoadCallback(function() { 

    setInterval(get_job_info, update_interval);
    // get_job_info();
    // get_aggregate_info();
});

function get_job_info() {
   jQuery.getJSON(domain_str + "names.json", process_job_info);
}


var job_info_data = [];

function process_job_info(data) {
    $("#program_name").text(data.program_name);
    $("#nprocs").text(data.nprocs + " processes");
    $("#current_time").text((data.time) + " seconds");

    // Render all the current values
    var container = $("#summary");
    var sorted_metrics = data.metrics.sort(function(a,b) { 
        return a.id - b.id; 
    });
    // Build an array of divs one for each metric with the name and value
    jQuery.each(sorted_metrics, function(i, metric) {
        var id = metric.id;
        var name = metric.name;
        var value = metric.rate;
        // if no job info has been created then create one as well as
        // the div to contain the display items
        if(job_info_data[id] == undefined) {
            // add a div the container
            var gauge_div_name = id + "_info_gauge";
            var str = 
                "<div class=\"metric_summary\" id=\"" + gauge_div_name  + "\">" +
                "<div class=\"name\">"  + name  + "</div>" +
                "<div class=\"value\">" + value + "</div>" +
                "<div class=\"gauge\"></div>" +
                "</div>";
            container.append(str);
            var div = container.children("#" + gauge_div_name);
            var gauge = new google.visualization.Gauge($(div).children(".gauge")[0]);
            job_info_data[id] = {
                div: div,
                gauge: gauge,
                options: {
                    width: 400, height: 120,
                    min: value, max: value + 1.0E-5},
                data:  google.visualization.arrayToDataTable([
                    ['Label', 'Value'], [name, value] ])
            };
        }
        // Get the job info
        var info = job_info_data[id];
        info.options.max = Math.max(info.options.max, value);
        info.options.min = Math.min(info.options.min, value);
        info.data.setCell(0,1, value);
        info.gauge.draw(info.data, info.options);
        info.div.children(".value").text(value);
    });
    // Get the job info again
}


function get_aggregate_info() {
    jQuery.getJSON(domain_str + "aggregate_all.json", process_aggregate_info);
}

function process_aggregate_info(data) {
    console.log(data);
    // Get the job info again
    setTimeout(process_aggregate_info, update_interval);
}



function render_applys(apply) {
    
    var data = google.visualization.arrayToDataTable([
        ['Label', 'Value'],
        ['Memory', 80],
        ['CPU', 55],
        ['Network', 68]
    ]);
    
}


// function print_info(data) {
//     console.log("printinfo");
//     console.log(data);
// }



/**
function drawChart() {
    var data = google.visualization.arrayToDataTable([
        ['Label', 'Value'],
        ['Memory', 80],
        ['CPU', 55],
        ['Network', 68]
    ]);
    
    var options = {
        width: 400, height: 120,
        redFrom: 90, redTo: 100,
        yellowFrom:75, yellowTo: 90,
        minorTicks: 5
    };
    
    var chart = new google.visualization.Gauge(document.getElementById('chart_div'));
    chart.draw(data, options);
}
*/

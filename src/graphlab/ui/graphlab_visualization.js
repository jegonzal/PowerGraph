google.load("jquery", "1.5");
google.load("jqueryui", "1.7.2");
google.load("visualization", "1", 
            {"packages":["corechart", "table", "gauge"]});


var domain_str = "http://localhost:8090/"
var update_interval = 1000;


// Start the rendering of the UI
google.setOnLoadCallback(function() { 
    initiate_job_info(); 
    initiate_aggregate_info();
});

function initiate_job_info() {
    var jqxhr = jQuery.getJSON(domain_str + "names.json", process_job_info)
        .error(function() { console.log("Unable to access " + domain_str + " will try again.");})
        .complete(function() {
            setTimeout(initiate_job_info, update_interval);
        });
}

function initiate_aggregate_info() {
    var jqxhr =    jQuery.getJSON(domain_str + "metrics_aggregate.json?rate=1", process_aggregate_info)
        .error(function() { console.log("Unable to access " + domain_str + " will try again.");})
        .complete(function() {
            setTimeout(initiate_aggregate_info, update_interval);
        });
}



var job_info_data = [];
function process_job_info(data) {
    console.log("Processing job info.");
    $("#program_name").text(data.program_name);
    // $("#nprocs").text(data.nprocs + " processes");
    $("#current_time").text((data.time) + " seconds");

    // Render all the current values
    var container = $("#gauges");
    var sorted_metrics = data.metrics.sort(function(a,b) { 
        return a.id - b.id; 
    });
    // Build an array of divs one for each metric with the name and value
    jQuery.each(sorted_metrics, function(i, metric) {
        var id = metric.id;
        var name = metric.name;
        var value = metric.rate_val;
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
}




var aggregate_charts = []

function process_aggregate_info(data) {
    console.log("Processing aggregate info.");

   // Render all the current values
    var container = $("#aggregate");
    var sorted_data = data.sort(function(a,b) { 
        return a.id - b.id; 
    });
    // Build an array of divs one for each metric with the name and value
    jQuery.each(sorted_data, function(i, metric) {
        var id = metric.id;
        var name = metric.name;
      
        if(aggregate_charts[id] == undefined) {
            // add a div the container
            var div_name = id + "_aggregate_chart";
            var str = 
                "<div class=\"aggregate\" id=\"" + div_name  + "\">" +
                "<div class=\"name\">"  + name  + "</div>" +
                "<div class=\"chart\"></div>" +
                "</div>";
            container.append(str);
            var div = $(container.children("#" +  div_name)).children(".chart")[0]; 
            aggregate_charts[id] = {
                div: div,
                options: { title: name, hAxis: {title: 'Time (seconds)',  
                                                titleTextStyle: {color: 'red'}}},
                chart: new google.visualization.AreaChart(div),
            }
        }
        // Update the chart
        var chart_info = aggregate_charts[id];
        chart_info.data = google.visualization.arrayToDataTable(
            [["Time", "Value"]].concat(metric.record));
        chart_info.chart.draw(chart_info.data, chart_info.options);
    });

}



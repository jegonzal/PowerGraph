google.load("jquery", "1.4.2");
google.load("jqueryui", "1.7.2");
google.load("visualization", "1");


var domain_str = "http://localhost:8090/wordclouds"
var update_interval = 1000;

var term_clouds = [];

// Start the rendering of the UI
google.setOnLoadCallback(function() { 
    setInterval(get_top_words, update_interval);
});

function get_top_words() {
    jQuery.getJSON(domain_str, process_top_words);
}




function process_top_words(data) {
    // Load summary info
    $("#ntopics").text(data.ntopics);
    $("#nwords").text(data.nwords);
    $("#ndocs").text(data.ndocs);
    $("#ntokens").text(data.ntokens);
    $("#alpha").text(data.alpha);
    $("#beta").text(data.beta);   

    // Render all the current values
    var container = $("#word_clouds");
    
    jQuery.each(data.values, function(i, term_count_table) {        
        if(term_clouds[i] == undefined) {
            var div_name = "term_cloud_" + i;            
            container.append(
                "<div class=\"cloud\" id=\"" + div_name  + "\"></div>");
            var div = container.children("#" + div_name);            
            var cloud = new TermCloud(div[0]);
            term_clouds[i] = { div: div, cloud: cloud };
        }
        var labels = [["String", "Value"]];        
        var table_data =  labels.concat(term_count_table);
        var table = google.visualization.arrayToDataTable(table_data);
        table.addColumn("string", "URL");
        //        console.log(table);
        term_clouds[i].cloud.draw(table, null );
    });
    // Get the job info again
} // end of process top words



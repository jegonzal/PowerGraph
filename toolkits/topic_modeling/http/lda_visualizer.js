google.load("jquery", "1.5");
google.load("jqueryui", "1.7.2");
google.load("visualization", "1");


var domain_str = "http://localhost:8090";
var page_str = "/wordclouds";
var update_interval = 2000;

function update_domain(form) {
    domain_str = form.inputbox.value;
    get_top_words();
}

var term_clouds = [];

// Start the rendering of the UI
google.setOnLoadCallback(function() { 
    get_top_words();
});


function get_top_words() {
    jQuery.getJSON(domain_str + page_str, process_top_words).error(function() { 
            console.log("Unable to access " + domain_str + " will try again.");
        }).complete(function() {
            setTimeout(get_top_words, update_interval);
        });
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



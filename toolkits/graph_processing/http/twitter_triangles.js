google.load("jquery", "1.5");
google.load("jqueryui", "1.7.2");
//google.load("visualization", "1");


var domain_str = "http://localhost:8090";
var domain_str = "";
var page_str =  "top_vertices_through.json";
var update_interval = 10000;
// jsonp callback required
var twitter_addr = "http://api.twitter.com/1/users/show.json"

function update_domain(form) {
    domain_str = form.inputbox.value;
    get_top_users();
}


// Start the rendering of the UI
google.setOnLoadCallback(function() { 
    get_top_users();
});


function get_top_users() {
    jQuery.getJSON(domain_str + page_str, process_top_users).error(function() { 
        console.log("Unable to access " + domain_str + " will try again.");
    });
    // .complete(function() {
    //             setTimeout(get_top_users, update_interval);
    //         });
}



var users_remaining = 0;
var top_user_profiles = [];

function process_top_users(data) {
    var user_ids = data.user_ids;
    users_remaining = user_ids.length;
    top_user_profiles = [];
    // launch jsonp requests for user data
    jQuery.each(user_ids, function(rank, user_id) {
        var query_str = twitter_addr + "?user_id=" + user_id; 
        jQuery.getJSON(query_str, function(data) { add_user(rank, data); }).error(function() { 
            console.log("Unable to access " + query_str + " will try again.");
        });
    });
}


function add_user(rank, data) {
    top_user_profiles[rank] = data;
    users_remaining--;
    if(users_remaining == 0) {
        render_page();
    }
}


function render_page() {
    var container = $("#top_degree");
    jQuery.each(top_user_profiles, function(rank, profile) {
        var div_name = profile.id_str;
        var div_contents = 
            "<div class=\"user\" id=\"" + profile.id_str + "\">" +
            "<img class=\"user_image\" src=\"" + profile.profile_image_url + "\" / >" +
            "<div class=\"name\">" +
            "  <a href=\"http://twitter.com/#!/" + profile.screen_name + "\">" + 
            "     " + profile.name + 
            "  </a>" +
            "</div>" +
            "</div>";
        container.append(div_contents);
    });
    

}


// function process_top_users(data) {
//     // Load summary info
//     $("#ntopics").text(data.ntopics);
//     $("#nwords").text(data.nwords);
//     $("#ndocs").text(data.ndocs);
//     $("#ntokens").text(data.ntokens);
//     $("#alpha").text(data.alpha);
//     $("#beta").text(data.beta);   

//     // Render all the current values
//     var container = $("#word_clouds");
    
//     jQuery.each(data.values, function(i, term_count_table) {        
//         if(term_clouds[i] == undefined) {
//             var div_name = "term_cloud_" + i;            
//             container.append(
//                 "<div class=\"cloud\" id=\"" + div_name  + "\"></div>");
//             var div = container.children("#" + div_name);            
//             var cloud = new TermCloud(div[0]);
//             term_clouds[i] = { div: div, cloud: cloud };
//         }
//         var labels = [["String", "Value"]];        
//         var table_data =  labels.concat(term_count_table);
//         var table = google.visualization.arrayToDataTable(table_data);
//         table.addColumn("string", "URL");
//         //        console.log(table);
//         term_clouds[i].cloud.draw(table, null );
//     });
//     // Get the job info again
// } // end of process top words



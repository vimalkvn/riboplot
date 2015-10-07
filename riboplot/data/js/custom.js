$(document).ready(function()
{
    $("#counts-table").tablesorter({sortList: [[1]]});
        //assign the sortStart event
    $("#counts-table").bind("sortStart", function() {
        $("#overlay").show();
    }).bind("sortEnd",function() {
        $("#overlay").hide();
    });
}
);

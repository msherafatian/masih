#' Include JavaScript for clipboard functionality
#'
#' @return JavaScript code as a string
#' @noRd
include_clipboard_js <- function() {
  tags$script(HTML("
    Shiny.addCustomMessageHandler('copyToClipboard', function(message) {
      navigator.clipboard.writeText(message.text).then(function() {
        Shiny.setInputValue(message.id, 'copied', {priority: 'event'});
      });
    });
  "))
}
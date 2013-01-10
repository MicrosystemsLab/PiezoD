<?php
//Function to output human readable file size
function format_size($size, $round = 0) {
    //Size must be bytes!
    $sizes = array('B', 'kB', 'MB', 'GB', 'TB', 'PB', 'EB', 'ZB', 'YB');
    for ($i=0; $size > 1024 && isset($sizes[$i+1]); $i++) $size /= 1024;
    return round($size,$round).$sizes[$i];
}

// Function to read the log file, and return an array as (filename => downloads)
function read_log()
{
	global $path;
	// Declare Array for holding data read from log file
	$name = array(); // array for file name
	$count = array(); // array for file count
	
	$file = @file("$path/Includes/log");
	if(empty($file))
	{
		return null;
	}
		
	// Read the entire contents of the log file into the arrays 
	$file = fopen("$path/Includes/log","r");
	while ($data = fscanf($file,"%[ -~]\t%d\n")) 
	{
		list ($temp1, $temp2) = $data;	
		array_push($name,$temp1);
		array_push($count,$temp2);
	}
	fclose($file);
	// $file_list contains data read from the log file as an array (filename => count)
	$file_list=@array_combine($name,$count); 
	ksort($file_list); // Sorting it in alphabetical order of key

	return $file_list;
}

// Reading the current directory to get a list of files and folders recursively
function scan_directory_recursively($directory, $filter=FALSE)
{
	// if the path has a slash at the end we remove it here
	if(substr($directory,-1) == '/')
	{
		$directory = substr($directory,0,-1);
	}

	// if the path is not valid or is not a directory ...
	if(!file_exists($directory) || !is_dir($directory))
	{
		// ... we return false and exit the function
		return FALSE;

	// ... else if the path is readable
	}elseif(is_readable($directory))
	{
		// we open the directory
		$directory_list = opendir($directory);

		// and scan through the items inside
		while (FALSE !== ($file = readdir($directory_list)))
		{
			// if the filepointer is not the current directory
			// or the parent directory
			if($file != '.' && $file != '..')
			{
				// we build the new path to scan
				$path = $directory.'/'.$file;

				// if the path is readable
				if(is_readable($path))
				{
					// we split the new path by directories
					$subdirectories = explode('/',$path);

					// if the new path is a directory
					if(is_dir($path))
					{
						// add the directory details to the file list
						$directory_tree[] = array(
							'path'    => $path,
							'name'    => end($subdirectories),
							'kind'    => 'directory',

							// we scan the new path by calling this function
							'content' => scan_directory_recursively($path, $filter));

					// if the new path is a file
					}elseif(is_file($path))
					{
						// get the file extension by taking everything after the last dot
						$extension = end(explode('.',end($subdirectories)));

						// if there is no filter set or the filter is set and matches
						if($filter === FALSE || $filter == $extension)
						{
							// add the file details to the file list
							$directory_tree[] = array(
								'path'      => $path,
								'name'      => end($subdirectories),
								'extension' => $extension,
								'size'      => filesize($path),
								'kind'      => 'file'
							);
						}
					}
				}
			}
		}
		// close the directory
		closedir($directory_list); 

		// return file list
		return @$directory_tree;

	// if the path is not readable ...
	}else{
		// ... we return false
		return FALSE;	
	}
}

// Function to recursively parse through the directory list and output html
function create_nested_list($items) 
{
	global $ignore;
	global $log;
	static $flag=0;
   
	if ($flag == 0) {
		$out = '<ul id="browser" class="filetree treeview-famfamfam">';
		$flag = 1;
	}
	else {
	   $out = '<ul>';
	}
	foreach ($items as $item) {
		if (in_array($item['name'],$ignore)) {
			continue;
		}	
		$out .= '<li>';
		if ($item['kind'] == 'directory') {
			$out .= '<span class="folder"><a><strong>' . $item['name'] . '</strong><em>-</em>-</a></span>';
			if ($item['content']) {
				$out .= create_nested_list($item['content']);
			}
		} 
		else {
			$out .= '<span class="' . strtolower($item['extension']) .'"><a href="download.php?fname='. $item['path'] .'"><strong>' .$item['name']. '</strong><em>' . format_size($item['size']) .'</em>';
			if (@array_key_exists($item['path'],$log)) {				
				$out .= '[' . $log[$item['path']] . ']';
			}
			else {
				$out .= '[0]';
			}
			$out .= '</a></span>';
		}
		$out .= "</li>";
   }
   
   $out .= "</ul>";
   return $out;
}

?>